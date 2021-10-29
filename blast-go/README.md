# Get GO annotation for non-model plant via blast

## Requirements
python3  
blast  
Good network for downloading  

## Usage
Configure the items in ./run_pipline.sh, while it is not necessary
```bash
bash run_blast_go.sh <input protein seq>.fasta
```
All required annotation database will be downloaded automatically if they cannot be detected
