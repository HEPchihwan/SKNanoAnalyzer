#!/usr/bin/env python3

import sys
import os

print("Python version:", sys.version)
print("Python executable:", sys.executable)
print("PYTHONPATH:", os.environ.get('PYTHONPATH', 'Not set'))

try:
    import ROOT
    print("✓ ROOT imported successfully")
    print("ROOT version:", ROOT.gROOT.GetVersion())
    
    # Test file access
    test_file = "/gv0/Users/achihwan/SKNanoRunlog/out/TTbar_test/2022/Muon_C.root"
    if os.path.exists(test_file):
        print(f"✓ Test file exists: {test_file}")
        
        try:
            f = ROOT.TFile.Open(test_file)
            if f and not f.IsZombie():
                print("✓ Can open ROOT file")
                f.Close()
            else:
                print("✗ Cannot open ROOT file")
        except Exception as e:
            print(f"✗ Error opening ROOT file: {e}")
    else:
        print(f"✗ Test file does not exist: {test_file}")
        
except ImportError as e:
    print(f"✗ Cannot import ROOT: {e}")
    print("Available modules in site-packages:")
    import site
    for path in site.getsitepackages():
        if os.path.exists(path):
            print(f"  {path}: {os.listdir(path)[:10]}")  # First 10 entries

try:
    import cmsstyle
    print("✓ cmsstyle imported successfully")
except ImportError as e:
    print(f"✗ Cannot import cmsstyle: {e}")

print("\nEnvironment variables:")
for key in ['ROOT_CONFIG_PATH', 'ROOTSYS', 'LD_LIBRARY_PATH']:
    print(f"  {key}: {os.environ.get(key, 'Not set')}")