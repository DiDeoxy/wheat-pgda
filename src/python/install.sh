# sneeds to be done all at once since source will exit when the script finishes
python3 -m venv src/python/ALIGN
source src/python/ALIGN/bin/activate
python -m pip install --upgrade pip
python -m pip install biopython pandas
deactivate