#!/bin/bash
# Quick Fix Script for TTN Agent Dependencies
# Fixes huggingface-hub version conflict

echo "=========================================="
echo "TTN Agent - Quick Dependency Fix"
echo "=========================================="
echo ""

# Activate virtual environment


echo ""
echo "📦 Fixing dependency conflicts..."
echo ""

# Fix huggingface-hub version conflict
echo "1. Downgrading huggingface-hub to compatible version..."
pip install "huggingface-hub>=0.30.0,<1.0" --force-reinstall

# Ensure transformers is up to date
echo ""
echo "2. Updating transformers..."
pip install transformers --upgrade

# Install vLLM if not present (optional)
echo ""
read -p "Do you want to install vLLM for local LLM support? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "3. Installing vLLM..."
    pip install vllm
else
    echo "3. Skipping vLLM installation"
fi

echo ""
echo "=========================================="
echo "✅ Dependencies fixed!"
echo "=========================================="
echo ""
echo "📝 Testing installation..."
python3 -c "
import sys
print('✓ Python imports check:')
try:
    import transformers
    print(f'  ✓ transformers {transformers.__version__}')
except Exception as e:
    print(f'  ✗ transformers: {e}')

try:
    from huggingface_hub import __version__
    print(f'  ✓ huggingface-hub {__version__}')
except Exception as e:
    print(f'  ✗ huggingface-hub: {e}')

try:
    import torch
    print(f'  ✓ torch {torch.__version__}')
except Exception as e:
    print(f'  ✗ torch: {e}')

try:
    import vllm
    print(f'  ✓ vLLM {vllm.__version__}')
except:
    print('  ℹ️  vLLM not installed (optional)')
"

echo ""
echo "=========================================="
echo "🎯 Next steps:"
echo "=========================================="
echo ""
echo "Test your installation:"
echo "  python main.py 2-178527548-C-T"
echo ""
echo "If you still see errors, try:"
echo "  pip install transformers --upgrade --force-reinstall"
echo ""
