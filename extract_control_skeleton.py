import argparse
import os


parser = argparse.ArgumentParser()
parser.add_argument(
    '--trait',
    type=str,
    required=True,
    help="Trait to be processed."
)
parser.add_argument(
    '--threshold',
    type=int,
    required=True,
    help="Threshold chosen for skeleton extraction."
)

args = parser.parse_args()

current_directory = os.path.dirname(__file__)

assert os.path.exists(os.path.join(current_directory, f'./output/{args.trait}/delete_control_score.txt')), print(
    'Must first run '
    'score_edges.py to obtain '
    'delete_control_score.txt')

tmp_folder = os.path.join(current_directory, f'./tmp_output/{args.trait}')
print('Extracting Control Skeleton...')
coreNodes = list(pd.read_table(os.path.join(tmp_folder, "P1_coreNodes.txt"), header=None)[0])
P1Nodes = list(pd.read_table(os.path.join(tmp_folder, "P1_neighborNodes.txt"), header=None)[0])
P2Nodes = list(pd.read_table(os.path.join(tmp_folder, "P2_neighborNodes.txt"), header=None)[0])
P3Nodes = list(pd.read_table(os.path.join(tmp_folder, "P3_neighborNodes.txt"), header=None)[0])
P4Nodes = list(pd.read_table(os.path.join(tmp_folder, "P4_neighborNodes.txt"), header=None)[0])

types = {}

edgeType_path = os.path.join(tmp_folder, 'control_network.txt')
edgeType_file = open(edgeType_path, 'r')
for line in edgeType_file:
    pas = line.strip().split()
    types.update({(pas[0], pas[2]): pas[1]})
edgeType_file.close()


controlDelete = []
delete_path = os.path.join(tmp_folder, 'delete_control_score.txt')
delete_file = open(delete_path, 'r')
for line in delete_file:
    pas = line.strip().split()
    controlDelete.append(float(pas[0]))
delete_file.close()

edgeScore = {}
edge_path = os.path.join(tmp_folder, 'edge_score.txt')
edge_file = open(edge_path, 'r')
for line in edge_file:
    pas = line.strip().split()
    edgeScore.update({(pas[0], pas[1]): (float(pas[2]), float(pas[3]))})
edge_file.close()

ascend_score = sorted(edgeScore.items(), key=lambda x: x[1][0], reverse=False)

output_path = os.path.join(current_directory, f'./output/{args.trait}/edges_delete_info.txt')
file = open(output_path, 'w')
for i in range(len(ascend_score)):
    strr = str(ascend_score[i][0][0]) + '\t' + str(types[ascend_score[i][0]]) + '\t' + str(
        ascend_score[i][0][1]) + '\t' + str(ascend_score[i][1][0]) + '\t' + str(controlDelete[i])
    file.write(strr)
    file.write('\n')
file.close()

controlSkeleton = []
path = os.path.join(current_directory, f'./output/{args.trait}/edges_delete_info.txt')
file = open(path, 'r')
for line in file:
    pas = line.strip().split()
    controlSkeleton.append((pas[0], pas[1], pas[2]))
file.close()

output_path = os.path.join(current_directory, f'./output/{args.trait}/control_skeleton.txt')
file = open(output_path, 'w')
for edge in controlSkeleton[15603:]:
    strr = str(edge[0]) + '\t' + str(edge[1]) + '\t' + str(edge[2])
    file.write(strr)
    file.write('\n')
file.close()

print('Done.')
print(f'Control Skeleton saved to ./output/{args.trait}/control_skeleton.txt')