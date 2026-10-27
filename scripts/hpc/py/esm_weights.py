"""Pre-fetch a Hugging Face model so GPU jobs can run on nodes without internet."""

import argparse

from huggingface_hub import snapshot_download


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--model", default="facebook/esm2_t33_650M_UR50D")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    path = snapshot_download(args.model, local_dir=args.out_dir,
                             allow_patterns=["*.json", "*.txt", "*.safetensors"])
    print(f"{args.model} -> {path}")


if __name__ == "__main__":
    main()
