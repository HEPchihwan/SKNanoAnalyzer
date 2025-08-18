# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

SKNanoAnalyzer is a C++ physics analysis framework for CMS (Compact Muon Solenoid) experiment data analysis. It processes NanoAOD files for Run 2 (2016-2018) and Run 3 (2022-2023) data periods, providing tools for dilepton studies, ttbar reconstruction, and systematic uncertainty evaluations.

## Build Commands

```bash
# Initial setup (required once per session)
source setup.sh

# Build entire project
./scripts/build.sh

# Incremental rebuild
./scripts/rebuild.sh
```

## Running Analysis

```bash
# Submit analysis jobs
SKNano.py -a AnalyzerName -i 'SamplePattern*' -e 2022 -n 10

# Example: Run DiLepton analyzer on DYJets samples
SKNano.py -a DiLepton -i 'DYJets*' -e 2022 -n 50 --reduction 10

# Skimming mode for data reduction
SKNano.py -a Skim_AnalyzerName -i SampleName -e 2022 --skimming_mode
```

## Architecture Overview

### Core Components
- **DataFormats/**: Physics object classes (Electron, Muon, Jet, etc.) with analysis methods
- **AnalyzerTools/**: Analysis utilities, corrections, systematic helpers, ML integration  
- **Analyzers/**: Specific analysis implementations inheriting from AnalyzerCore
- **External dependencies**: KinematicFitter, RoccoR, jsonpog-integration, LHAPDF, libtorch

### Key Classes
- **AnalyzerCore**: Base class providing common functionality for all analyzers
- **SKNanoLoader**: Handles NanoAOD file loading and branch management
- **SystematicHelper**: Manages systematic uncertainty variations
- **MyCorrection**: Centralized correction and scale factor application

### Data Organization
- Era-based structure: `data/{era}/` contains corrections, efficiencies, and sample metadata
- Sample information in `data/{era}/Sample/CommonSampleInfo.json`
- Systematic configurations in YAML files (DiLeptonSystematic.yaml, VcbSystematic.yaml)

## Development Patterns

### Creating New Analyzers
1. Create header in `Analyzers/include/` inheriting from `AnalyzerCore`
2. Implement in `Analyzers/src/` with required methods:
   - `initializeAnalyzer()`
   - `executeEvent()`
   - Virtual destructor
3. Add to `Analyzers/include/AnalyzersLinkDef.hpp`
4. Use enum classes for constants instead of magic numbers/strings

### Systematic Uncertainties
- Use `SystematicHelper` class for systematic variations
- Implement variations using `variation` enum (nom, up, down)
- Apply corrections through `MyCorrection` class methods

### Physics Object Usage
```cpp
// Example pattern for lepton selection
RVec<Lepton *> leptons;
TString MuonID = "Tight"; // Example ID
for (auto &mu : *muons) {
    if (mu.Pass(MuonID)) {
        leptons.push_back(&mu);
    }
}
```

## Dependencies and Environment

### Required Packages (managed via conda/mamba)
- ROOT 6.32+
- correctionlib (CMS corrections)
- onnxruntime-cpp (ML inference)
- PyTorch/libtorch (deep learning)
- LHAPDF (auto-installed via scripts/install\_lhapdf.sh)

### Configuration
- Personal config: `config/config.$USER` (system type, package manager, notifications)
- Sample lists: `SampleLists/{era}/` for different analysis categories

## Job Management

### HTCondor Integration
- Jobs submitted via HTCondor with DAG workflows
- Configurable memory (default 2GB) and CPU resources
- Automatic output file management and merging
- Telegram bot notifications (optional)

### Sample Selection
- Regex patterns for flexible sample matching
- Era-specific sample organization (2016preVFP, 2017, 2018, 2022, 2022EE, 2023, etc.)
- Cross-section and luminosity handled automatically

## Multi-Era Support

The framework supports unified analysis across data-taking periods:
- **Run 2**: 2016preVFP, 2016postVFP, 2017, 2018
- **Run 3**: 2022, 2022EE, 2023, 2023BPix

Era-specific corrections and configurations are automatically loaded based on the `-e` flag in job submission.

## Environment Setup and Troubleshooting

### Common Issues and Solutions

#### ROOT and Library Conflicts
- **GLIBCXX Version Conflicts**: If ROOT commands (like `hadd`) fail with GLIBCXX version errors:
  - The issue occurs when conda ROOT requires newer GLIBCXX than system libstdc++
  - Solution: setup.sh automatically prioritizes conda environment libraries in LD_LIBRARY_PATH
  - A libstdc++.so.6 symlink is created in conda environment if needed

#### Missing Dependencies
- **HTCondor Python Module**: Install with `micromamba install -y htcondor`
  - Provides python-htcondor, htcondor-cli, htcondor-utils packages
  - Required for job submission and management

- **Boost Libraries**: Install complete Boost package:
  ```bash
  micromamba install -c conda-forge boost-cpp libboost-devel
  ```
  - Provides headers in `$CONDA_PREFIX/include/boost/`
  - CMakeLists.txt automatically detects ROOT's built-in Boost or conda Boost

#### CMake Configuration Issues
- **correctionlib Not Found**: CMakeLists.txt sets correctionlib path explicitly:
  ```cmake
  set(correctionlib_DIR "/path/to/conda/lib/python3.9/site-packages/correctionlib/cmake")
  ```
  - Path is automatically configured based on conda environment
  - Use `correction config --cmake` to get the correct path

- **CMAKE_PREFIX_PATH**: Automatically includes conda environment:
  ```cmake
  set(CMAKE_PREFIX_PATH "$CONDA_PREFIX;${CMAKE_PREFIX_PATH}")
  ```

- **yaml-cpp ABI Compatibility Issues**: If build fails with yaml-cpp linking errors:
  ```bash
  # Install conda yaml-cpp to avoid ABI conflicts
  micromamba install -y yaml-cpp
  ```
  - CMakeLists.txt automatically prefers system yaml-cpp over FetchContent
  - AnalyzerTools links to `yaml-cpp` (not `yaml-cpp::yaml-cpp`) for compatibility
  - Clean build directory if switching from FetchContent to system version

### Environment Variables Set by setup.sh
- **CONDA_PREFIX**: Explicitly set to conda environment path for CMake
- **LD_LIBRARY_PATH**: Prioritizes conda libraries over system libraries
- **CORRECTION_INCLUDE_DIR/CORRECTION_LIB_DIR**: For correctionlib integration
- **ONNXRUNTIME_INCLUDE_DIR/ONNXRUNTIME_LIB_DIR**: For ML inference
- **LHAPDF_***: Paths for PDF sets and libraries

### Build System Enhancements
- **Automatic Boost Detection**: CMakeLists.txt checks multiple sources:
  1. ROOT's built-in Boost (preferred)
  2. Conda environment Boost
  3. System Boost (fallback)
- **Library Priority**: conda libraries take precedence to avoid version conflicts
- **Multi-threading**: Uses `make -j6` for parallel compilation

### Verification Commands
```bash
# Test ROOT functionality
source setup.sh && hadd --help

# Verify HTCondor Python module
source setup.sh && python -c "import htcondor; print('HTCondor version:', htcondor.version())"

# Check cmake dependencies
source setup.sh && correction config --cmake

# Test build system
source setup.sh && ./scripts/build.sh
```

### Notes for Future Maintenance
- Always run `source setup.sh` before building or running analyses
- The setup script handles conda environment activation and library path configuration
- CMake configuration is robust against missing environment variables
- All paths are automatically detected based on the active conda environment
