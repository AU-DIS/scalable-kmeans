alias python3="$(which python3)"

filename = $1

PYCMD=$(cat <<EOF
import h5py
import numpy as np
import pandas as pd
import itertools as it;

print('$1')
f = h5py.File('$1', 'r')
print(list(f.keys()))
dataset = f['data']

arr = np.array(dataset)

df = pd.DataFrame(data=[(*axes, v.item()) for axes, v in zip(it.product(*[range(i) for i in arr.shape]), np.nditer(arr))], columns=tuple('xyzkv'))

#print(df)
pd.DataFrame(df).to_csv('$1.txt', header=False, index=False, sep=" ", line_terminator=" \n")

EOF
)

TEMP_SCRIPT=$(mktemp)
echo "$PYCMD" > "$TEMP_SCRIPT"
python3 "$TEMP_SCRIPT"