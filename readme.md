# CSon Pipeline

![workflow](workflow.bmp)

## Overview

This pipeline identifies control skeletons in omnigenic networks through a 3-stage computational process. The workflow calculates edge criticality, determines optimal thresholds, and extracts core control structures for disease trait analysis.

## File Descriptions

### Core Modules
| File                          | Functionality                                                |
| ----------------------------- | ------------------------------------------------------------ |
| `score_edges.py`              | Calculates optimal control structures and edge criticality scores |
| `choose_threshold.py`         | Visualizes network control capacity curve and suggests thresholds |
| `extract_control_skeleton.py` | Extracts final control skeleton using determined threshold   |
| `utils/ms.py`                 | Maximum matching generation                                  |
| `utils/utils.py`              | Shared utility functions and I/O operations                  |

### Input Requirements (per trait)
```text
data/
  └─ [trait_name]/
      ├─ omnigenic network.txt
      ├─ genetic perturbation values.txt
      ├─ omnigenic neighborhood/
          ├─ core.txt
          └─ periphery.txt
```

## Execution Pipeline

### Step 1: Edge Scoring
```bash
python score_edges.py --trait [trait_name]
```
**Outputs:**
- `output/[trait_name]/delete_control_score.txt`
- `output/[trait_name]/edge_score.txt`

### Step 2: Threshold Determination
```bash
python choose_threshold.py --trait [trait_name]
```
**Output:** Interactive threshold selection plot

### Step 3: Skeleton Extraction
```bash
python extract_control_skeleton.py --trait [trait_name] --threshold [N]
```
**Final Output:**
- `output/[trait_name]/control_skeleton.txt`

## Custom Trait Integration
To analyze new traits:
1. Create trait directory under `data/`
2. Provide required input files as mentioned in 'Input Requirements'
3. Follow standard execution pipeline

