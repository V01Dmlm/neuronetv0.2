import torch
from transformers import AutoTokenizer, AutoModelForCausalLM
from peft import PeftModel
from pathlib import Path

def load_model(model_dir="D:/projects/neuronet-v0.1-main/neuronet-v0.1-main/gptneo-1.3B-lora-finetuned"):
    model_path = Path(model_dir)
    if not model_path.exists():
        raise ValueError(f"Model directory does not exist: {model_path}")

    # Base model (pretrained GPT-Neo 1.3B)
    base_model_name = "EleutherAI/gpt-neo-1.3B"
    print(f"Loading base model from {base_model_name}...")
    base_model = AutoModelForCausalLM.from_pretrained(base_model_name)

    # Tokenizer (from the fine-tuned folder if available)
    print(f"Loading tokenizer from {model_dir}...")
    tokenizer = AutoTokenizer.from_pretrained(model_dir, local_files_only=True)
    tokenizer.pad_token = tokenizer.eos_token  # Make sure padding token is defined

    # Load LoRA adapter
    print(f"Applying LoRA adapter from {model_dir}...")
    model = PeftModel.from_pretrained(base_model, model_dir, torch_dtype=torch.float16, local_files_only=True)

    model.eval()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.to(device)
    print(f"Model loaded on {device}")
    return model, tokenizer, device


def generate_summary(model, tokenizer, device, abstract_text, max_length=150):
    prompt = f"Summarize the following abstract:\n{abstract_text}\nSummary:"
    inputs = tokenizer(prompt, return_tensors="pt", truncation=True, max_length=512)
    input_ids = inputs.input_ids.to(device)
    attention_mask = inputs.attention_mask.to(device)

    with torch.no_grad():
        # Call generate on the underlying base_model to avoid PEFT generate error
        outputs = model.base_model.generate(
            input_ids=input_ids,
            attention_mask=attention_mask,
            max_length=max_length + input_ids.shape[1],
            do_sample=True,
            temperature=0.7,
            top_p=0.9,
            eos_token_id=tokenizer.eos_token_id,
            pad_token_id=tokenizer.eos_token_id,
            num_return_sequences=1
        )

    summary = tokenizer.decode(outputs[0][input_ids.shape[1]:], skip_special_tokens=True)
    return summary.strip()


if __name__ == "__main__":
    model, tokenizer, device = load_model()
    sample_abstract = (
        "Alzheimer's disease is a neurodegenerative disorder "
        "characterized by progressive memory loss and cognitive decline..."
    )
    print("Generating summary...")
    summary = generate_summary(model, tokenizer, device, sample_abstract)
    print("Summary:\n", summary)
