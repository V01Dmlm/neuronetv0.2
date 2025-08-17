import sys
import threading
from collections import deque

from PyQt6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel,
    QTextEdit, QLineEdit, QPushButton, QSpinBox, QDoubleSpinBox,
    QComboBox, QSizePolicy, QDialog, QScrollArea
)
from PyQt6.QtCore import Qt

from generate_pubmed_dataset import generate_dataset as scrape_main
from finetune import main as train_main
from generate_summary import load_model, generate_summary
from pubmed_fetcher import fetch_top_abstract 

# ===================== Dark Theme =====================
def apply_dark_theme(app):
    dark_stylesheet = """
        QWidget {
            background-color: #121212;
            color: #ffffff;
            font-size: 14px;
        }
        QLineEdit, QTextEdit, QSpinBox, QDoubleSpinBox, QComboBox {
            background-color: #1e1e1e;
            color: #ffffff;
            border: 1px solid #444;
        }
        QPushButton {
            background-color: #2e2e2e;
            color: #fff;
            border: 1px solid #555;
            padding: 6px;
        }
        QPushButton:hover {
            background-color: #444;
        }
        QLabel {
            color: #ccc;
        }
    """
    app.setStyleSheet(dark_stylesheet)

MAX_HISTORY = 6
chat_history = deque(maxlen=MAX_HISTORY)

# ===================== Settings Pop-up =====================
class SettingsDialog(QDialog):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("⚙️ Settings")
        self.setGeometry(200, 200, 500, 600)

        layout = QVBoxLayout()

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        inner_widget = QWidget()
        inner_layout = QVBoxLayout()

        # Dataset section
        dataset_layout = QVBoxLayout()
        self.query_input = QLineEdit("Alzheimer OR Cancer OR Diabetes")
        self.max_results_input = QSpinBox()
        self.max_results_input.setMaximum(10000)
        self.max_results_input.setValue(1000)
        self.prompt_input = QTextEdit()
        self.prompt_input.setFixedHeight(100)
        self.prompt_input.setPlainText("Summarize the following abstract:\n{abstract}\nAnswer:")
        dataset_btn = QPushButton("Download & Generate Dataset")
        dataset_btn.clicked.connect(self.run_scraper_thread)

        dataset_layout.addWidget(QLabel("📄 Generate Dataset"))
        dataset_layout.addWidget(QLabel("Query:"))
        dataset_layout.addWidget(self.query_input)
        dataset_layout.addWidget(QLabel("Max Results:"))
        dataset_layout.addWidget(self.max_results_input)
        dataset_layout.addWidget(QLabel("Prompt Template:"))
        dataset_layout.addWidget(self.prompt_input)
        dataset_layout.addWidget(dataset_btn)

        # Training section
        train_layout = QVBoxLayout()
        self.epochs_input = QSpinBox()
        self.epochs_input.setMinimum(1)
        self.epochs_input.setValue(5)
        self.lr_input = QDoubleSpinBox()
        self.lr_input.setDecimals(6)
        self.lr_input.setValue(3e-5)
        self.lr_input.setSingleStep(1e-5)
        train_btn = QPushButton("Train model")
        train_btn.clicked.connect(self.run_training_thread)

        train_layout.addWidget(QLabel("🛠️ Train Model"))
        train_layout.addWidget(QLabel("Epochs:"))
        train_layout.addWidget(self.epochs_input)
        train_layout.addWidget(QLabel("Learning Rate:"))
        train_layout.addWidget(self.lr_input)
        train_layout.addWidget(train_btn)

        # Combine sections
        inner_layout.addLayout(dataset_layout)
        inner_layout.addSpacing(20)
        inner_layout.addLayout(train_layout)
        inner_layout.addStretch()
        inner_widget.setLayout(inner_layout)
        scroll.setWidget(inner_widget)
        layout.addWidget(scroll)
        self.setLayout(layout)

    def run_scraper_thread(self):
        def scraper_task():
            query = self.query_input.text()
            max_results = self.max_results_input.value()
            prompt_template = self.prompt_input.toPlainText()
            output_file = "neuronet_data.jsonl"
            print(f"Scraping {query} ...")
            try:
                scrape_main(query=query, max_results=max_results,
                            output_file=output_file,
                            prompt_template=prompt_template)
                print(f"✅ Dataset saved to {output_file}")
            except Exception as e:
                print(f"❌ Error: {e}")
        threading.Thread(target=scraper_task).start()

    def run_training_thread(self):
        def train_task():
            print("🚀 Starting training...")
            try:
                train_main()
                print("✅ Training complete!")
            except Exception as e:
                print(f"❌ Error: {e}")
        threading.Thread(target=train_task).start()

# ===================== Main GUI =====================
class PubMedSummarizerGUI(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("🧠 NeuroNet - PubMed Research Assistant")
        self.setGeometry(100, 100, 900, 650)

        main_layout = QVBoxLayout()

        # ===================== Top Bar =====================
        top_bar = QHBoxLayout()
        self.mood_combo = QComboBox()
        self.mood_combo.addItems([
            "🧠 Researcher Mode",
            "🩺 Clinical Summary",
            "👨‍🏫 Professor Mode",
            "💬 Casual Explainer"
        ])
        self.mood_combo.setFixedWidth(200)
        top_bar.addWidget(self.mood_combo)
        top_bar.addStretch()  # push settings to right

        # Settings gear button
        self.settings_button = QPushButton("⚙️")
        self.settings_button.setFixedWidth(50)
        self.settings_button.clicked.connect(self.open_settings_popup)
        top_bar.addWidget(self.settings_button)
        main_layout.addLayout(top_bar)

        # ===================== Chat Display =====================
        self.chat_display = QTextEdit()
        self.chat_display.setReadOnly(True)
        self.chat_display.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        main_layout.addWidget(self.chat_display)

        # ===================== Input Row =====================
        input_layout = QHBoxLayout()
        self.user_input = QLineEdit()
        self.user_input.setPlaceholderText("Ask a medical research question...")
        self.user_input.returnPressed.connect(self.chatbot_answer)

        self.send_button = QPushButton("Send")
        self.send_button.setFixedWidth(80)
        self.send_button.clicked.connect(self.chatbot_answer)

        input_layout.addWidget(self.user_input)
        input_layout.addWidget(self.send_button)
        main_layout.addLayout(input_layout)

        self.setLayout(main_layout)

        # ===================== Welcome Message =====================
        self.chat_display.append("🤖 Welcome, researcher! How can I help you today?")

        # ===================== Style Templates =====================
        self.style_map = {
            "🧠 Researcher Mode": (
                "You are NeuroNet, a smart and curious medical research assistant.\n\n"
                "Abstract:\n{abstract}\n\n"
                "User: {user_question}\nNeuroNet:"
            ),
            "🩺 Clinical Summary": (
                "Summarize this clinical abstract clearly, with focus on results.\n\n"
                "Abstract:\n{abstract}\n\nSummary:"
            ),
            "👨‍🏫 Professor Mode": (
                "You're a professor. Explain this abstract to a curious med student.\n\n"
                "Abstract:\n{abstract}\n\n"
                "Student: {user_question}\nProfessor:"
            ),
            "💬 Casual Explainer": (
                "Explain the abstract in simple language.\n\n"
                "Abstract:\n{abstract}\n\n"
                "User: {user_question}\nYou:"
            )
        }

    # ===================== Open Settings Popup =====================
    def open_settings_popup(self):
        self.settings_dialog = SettingsDialog()
        self.settings_dialog.exec()

    # ===================== Chatbot =====================
    def chatbot_answer(self):
        user_msg = self.user_input.text().strip()
        if not user_msg:
            return

        selected_style_name = self.mood_combo.currentText()
        template = self.style_map[selected_style_name]

        self.chat_display.append(f"👤 You: {user_msg}")
        self.user_input.clear()

        def generate():
            try:
                self.chat_display.append("📡 Fetching related abstract from PubMed...")
                abstract = fetch_top_abstract(user_msg)
                if not abstract:
                    self.chat_display.append("❌ No abstract found.")
                    return

                full_prompt = template.replace("{abstract}", abstract)\
                                    .replace("{user_question}", user_msg)

                model, tokenizer, device = load_model()
                answer = generate_summary(model, tokenizer, device, full_prompt)
                self.chat_display.append(f"🤖 NeuroNet:\n{answer}")
                chat_history.append(f"User: {user_msg}")
                chat_history.append(f"NeuroNet: {answer}")
            except Exception as e:
                self.chat_display.append(f"[Error] {e}")

        threading.Thread(target=generate).start()

# ===================== App Entry =====================
if __name__ == "__main__":
    app = QApplication(sys.argv)
    apply_dark_theme(app)
    window = PubMedSummarizerGUI()
    window.show()
    sys.exit(app.exec())
