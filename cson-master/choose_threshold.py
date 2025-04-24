import argparse
import os
import numpy as np

from matplotlib import pyplot as plt
import matplotlib.cm as cm

from utils.utils import find_max_drop_index, RearKNormalize

parser = argparse.ArgumentParser()
parser.add_argument(
    '--trait',
    type=str,
    required=True,
    help="Trait to be processed."
)

args = parser.parse_args()

current_directory = os.path.dirname(__file__)

assert os.path.exists(os.path.join(current_directory, f'./output/{args.trait}/delete_control_score.txt')), print(
    'Must first run '
    'score_edges.py to obtain '
    'delete_control_score.txt')

controlDelete = []
score_path = os.path.join(current_directory, f'./output/{args.trait}/delete_control_score.txt')
score_file = open(score_path, 'r')
for line in score_file:
    pas = line.strip().split()
    controlDelete.append(float(pas[0]))
score_file.close()

threshold = find_max_drop_index(controlDelete)

print(f'Suggested threshold: {threshold}')
print('Plotting curve...')

plt.rcParams["font.family"] = "Times New Roman"
fig, ax = plt.subplots(figsize=(5, 5))

colormap = cm.YlGn
y = controlDelete
x = np.array([i for i in range(len(y))])
norm = RearKNormalize(vmin=np.min(x), vmax=np.max(x), rear_k=0.7)

for i in range(len(x) - 1):
    ax.plot(x[i:i+2], y[i:i+2], marker='s', markersize=1.5, linewidth=1, color=colormap(norm(x[i])))


ax.vlines(threshold, -1, 18.164, color='r', linestyle='--')

ax.set_ylim(-1, 44)
ax.set_xticks([threshold])
ax.set_xticklabels([threshold])
ax.set_yticklabels([])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

plt.title(f'{args.trait} threshold choosing')
plt.show()
