import argparse
import os
import networkx as nx
import pandas as pd

from tqdm import tqdm

from utils.utils import *
from utils.ms import enum_and_save_Ms

parser = argparse.ArgumentParser()
parser.add_argument(
    '--trait',
    type=str,
    required=True,
    help="Trait to be processed."
)

args = parser.parse_args()

current_directory = os.path.dirname(__file__)
assert os.path.exists(os.path.join(current_directory, f'./data/{args.trait}/')), print(
    'Pls first create trait folder in data')
tmp_folder = os.path.join(current_directory, f'./tmp_output/{args.trait}')

if not os.path.exists(tmp_folder):
    os.mkdir(tmp_folder)

with open(os.path.join(current_directory, f'./data/{args.trait}/omnigenic network.txt'), 'r') as input_file:
    with open(os.path.join(tmp_folder, f'before_binary.txt'), 'w') as output1, \
            open(os.path.join(tmp_folder, f'binary.txt'), 'w') as output2:
        for line in input_file.readlines():
            output1.write(line)
            gene_1, gene_2 = line.strip('\n').split('\t')
            output2.write(gene_1 + '_out\t' + gene_2 + '_in\n')

print("[1 / 10] Loading omnigenic network...")
allDG = nx.DiGraph()
network_path = os.path.join(tmp_folder, 'before_binary.txt')
f = open(network_path, 'r')
for line in f:
    pas = line.strip().split()
    allDG.add_node(pas[0])
    allDG.add_node(pas[1])
    allDG.add_edge(pas[0], pas[1])
f.close()

print("[2 / 10] Loading core nodes...")
core_nodes = list(
    pd.read_table(os.path.join(current_directory, f'./data/{args.trait}/omnigenic neighborhood/core.txt'), header=None)[
        0])

tempNodes = set(core_nodes)
neighbors, neighborEdges = getNeighborsFromNodes(allDG, tempNodes)

print('[3 / 10] Stratifying network...')

index = 1
while ((index == 1) or (len(neighbors - tempNodes) >= 600)):
    neighbors = neighbors - tempNodes
    coreEdges = list(allDG.subgraph(list(tempNodes)).edges)  # get core edges
    coreNodes = list(allDG.subgraph(list(tempNodes)).nodes)  # core nodes
    neighborNodes = neighbors - tempNodes  # neighbor of core
    subEdges = coreEdges + neighborEdges
    subEdges = set(subEdges)
    # delete backlinks
    deleteEdges = []
    for edge in subEdges:
        if edge[0] in coreNodes and edge[1] not in coreNodes:
            deleteEdges.append(edge)
    subEdges = list(set(subEdges) - set(deleteEdges))
    # save core nodes
    coreNodes_path = os.path.join(tmp_folder, f'P{index}_coreNodes.txt')
    file = open(coreNodes_path, 'w')
    strr = ''
    for node in coreNodes:
        strr = str(node)
        file.write(strr)
        file.write('\n')
    file.close()
    # save core neighbors
    neighborNodes_path = os.path.join(tmp_folder, f'P{index}_neighborNodes.txt')
    file = open(neighborNodes_path, 'w')
    strr = ''
    for node in neighborNodes:
        strr = str(node)
        file.write(strr)
        file.write('\n')
    file.close()
    # save subgraph as bipartite
    output_path1 = os.path.join(tmp_folder, f'P{index}_binary.txt')
    output_path2 = os.path.join(tmp_folder, f'P{index}_non_binary.txt')
    file1 = open(output_path1, 'w')
    file2 = open(output_path2, 'w')
    strr1 = ''
    strr2 = ''
    for edge in subEdges:
        strr1 = str(edge[0]) + '_out' + '\t' + str(edge[1]) + '_in'
        strr2 = str(edge[0]) + '\t' + str(edge[1])
        file1.write(strr1)
        file2.write(strr2)
        file1.write('\n')
        file2.write('\n')
    file1.close()
    file2.close()
    tempNodes = tempNodes | neighbors
    neighbors, neighborEdges = getNeighborsFromNodes(allDG, neighbors)

    index += 1

neighborNodesRemain = neighbors - tempNodes
lastNeighborNodesRemain = set()
neighborEdgesRemain = neighborEdges
while (len(neighborNodesRemain - lastNeighborNodesRemain) != 0):
    lastNeighborNodesRemain = neighborNodesRemain
    neighbors, neighborEdges = getNeighborsFromNodes(allDG, lastNeighborNodesRemain)
    neighbors = neighbors - tempNodes
    neighborNodesRemain = lastNeighborNodesRemain | neighbors
    neighborEdgesRemain = neighborEdgesRemain + neighborEdges

coreEdges = list(allDG.subgraph(list(tempNodes)).edges)
coreNodes = list(allDG.subgraph(list(tempNodes)).nodes)
neighborNodes = neighborNodesRemain - tempNodes
subEdges = coreEdges + neighborEdgesRemain
subEdges = set(subEdges)
deleteEdges = []
for edge in subEdges:
    if edge[0] in coreNodes and edge[1] not in coreNodes:
        deleteEdges.append(edge)
subEdges = list(set(subEdges) - set(deleteEdges))
# save cores
coreNodes_path = os.path.join(tmp_folder, f'P{index}_coreNodes.txt')
file = open(coreNodes_path, 'w')
strr = ''
for node in coreNodes:
    strr = str(node)
    file.write(strr)
    file.write('\n')
file.close()
# save core neighbors
neighborNodes_path = os.path.join(tmp_folder, f'P{index}_neighborNodes.txt')

file = open(neighborNodes_path, 'w')
strr = ''
for node in neighborNodes:
    strr = str(node)
    file.write(strr)
    file.write('\n')
file.close()
# saving subgraph as bipartite
output_path1 = os.path.join(tmp_folder, f'P{index}_binary.txt')
output_path2 = os.path.join(tmp_folder, f'P{index}_non_binary.txt')
file1 = open(output_path1, 'w')
file2 = open(output_path2, 'w')
strr1 = ''
strr2 = ''
for edge in subEdges:
    strr1 = str(edge[0]) + '_out' + '\t' + str(edge[1]) + '_in'
    strr2 = str(edge[0]) + '\t' + str(edge[1])
    file1.write(strr1)
    file2.write(strr2)
    file1.write('\n')
    file2.write('\n')
file1.close()
file2.close()

print('[4 / 10] Simplifying network...')

for s in ['non_binary', 'binary']:
    for i in range(1, index):
        C2_path = os.path.join(tmp_folder, f'P{i + 1}_{s}.txt')
        C2_f = open(C2_path, 'r')
        C2_edges = set()
        for lines in C2_f:
            pas = lines.strip().split()
            C2_edges.add((pas[0], pas[1]))
        C2_f.close()

        C1_path = os.path.join(tmp_folder, f'P{i}_{s}.txt')
        C1_f = open(C1_path, 'r')
        C1_edges = set()
        for lines in C1_f:
            pas = lines.strip().split()
            C1_edges.add((pas[0], pas[1]))
        C1_f.close()
        C2_edges_selected = C2_edges - C1_edges

        output_path = os.path.join(tmp_folder, f'P{i + 1}_{s}_selected.txt')
        file = open(output_path, 'w')
        strr = ''
        for edge in C2_edges_selected:
            strr = str(edge[0]) + '\t' + str(edge[1])
            file.write(strr)
            file.write('\n')
        file.close()

# comment this 'if' when running
# if not os.path.exists(os.path.join(tmp_folder, f'P1_maxmatching.txt')):
print('[5 / 10] Enumerating Ms... This is going to take a long time...')
enum_and_save_Ms(tmp_folder, index)

# Calculate domination centrality
# if not os.path.exists(os.path.join(tmp_folder, './ms')):
print('[6 / 10] Calculating domination centrality...')
os.mkdir(os.path.join(tmp_folder, './ms'))
for i in range(1, index + 1):
    os.mkdir(os.path.join(tmp_folder, f'./ms/P{i}_ms/'))

    path = os.path.join(tmp_folder, f'P{i}_maxmatching.txt')
    allList = MsToList(path)

    msList = allList
    allDG = nx.DiGraph()
    if i == 1:
        network_path = os.path.join(tmp_folder, f'P{i}_non_binary.txt')
    else:
        network_path = os.path.join(tmp_folder, f'P{i}_non_binary_selected.txt')

    allEdges = []
    f = open(network_path, 'r')
    for line in f:
        pas = line.strip().split()
        allDG.add_node(pas[0])
        allDG.add_node(pas[1])
        allDG.add_edge(pas[0], pas[1])
        allEdges.append((pas[0], pas[1]))
    f.close()

    # Calculate domination centrality
    # count of ms
    count = 0
    pbar = tqdm(msList, desc=f'Processing P{i}')
    for ms in pbar:
        # ms to nx.graph
        msDG = nx.DiGraph()
        msDG.add_edges_from(ms)
        # find cycle nodes of ms
        cycle = set()
        for c in sorted(nx.strongly_connected_components(msDG), key=len, reverse=True):
            if len(c) > 1:
                cycle = cycle | c
        # find stem nodes of ms
        cycleNodes = list(cycle)  # cycle node
        msNodes = msDG.nodes  # all node in ms
        stemNodes = list(msNodes - cycleNodes)  # stem node
        # Determine which vertices are in the STEM node set,
        # i.e., the degree of entry is 0, which cannot be used as the starting point of additional links
        # also, determine which are the tail points,
        # i.e., the out-degree is 0, and when the network is reversed for OS calculation, the node is the stem "vertex"
        stemTopNodes = []
        stemBottomNodes = []
        for node in stemNodes:
            if msDG.in_degree(node) == 0:
                stemTopNodes.append(node)
            if msDG.out_degree(node) == 0:
                stemBottomNodes.append(node)
        # Gets the set of non-vertex nodes in the stem nodes
        stemNotTopNodes = list(set(stemNodes) - set(stemTopNodes))
        # Obtain the set of nodes with non-ending points in the stem nodes
        stemNotBottomNodes = list(set(stemNodes) - set(stemBottomNodes))
        # Gets the largest set of non-matching nodes in the ms
        unMatchingNodes = list(allDG.nodes - msDG.nodes)
        # Gets the largest set of non-matching edges in the ms
        unMatchingEdges = list(set(allEdges) - set(ms))
        # select stem->cycle  cycle->stem from non-matching edges
        stemToCycleEdges = []
        cycleToStemEdges = []
        for edge in unMatchingEdges:
            if (edge[0] in stemNotTopNodes) and (edge[1] in cycleNodes):
                stemToCycleEdges.append(edge)
            if (edge[0] in cycleNodes) and (edge[1] in stemNotBottomNodes):
                cycleToStemEdges.append(edge)
        # add the additional edges of stem->cycle to the temporary CS maximum matching edge set
        tempCsMs = ms
        tempCsMs.extend(stemToCycleEdges)
        # used to calculate information about the control subspace of the STEM node
        tempCsMsDG = nx.DiGraph()
        tempCsMsDG.add_edges_from(tempCsMs)
        # add the additional edges of cycle->stem to the maximum matching edge set of the temporary OS
        tempOsMs = ms
        tempOsMs.extend(cycleToStemEdges)
        # used to calculate observational subspace-related information for STEM nodes
        tempOsMsDG = nx.DiGraph()
        tempOsMsDG.add_edges_from(tempOsMs)
        tempOsMsDG = tempOsMsDG.reverse(copy=True)  # reverse

        nodeInfo = []
        # every ms as a file
        for node in list(allDG.nodes):

            if node in unMatchingNodes:
                nodeInfo.append((node, -1, -1, -1, -1, -1, -1, -1))
            elif node in cycleNodes:  # node in ms, and is a cycle node
                # calculating observable subspace
                rankControl, controlNodes, controlEdges = getRank(msDG, node, False)
                # calculating observable subspace
                msDG = msDG.reverse(copy=True)  # reverse
                rankObserve, observeNodes, observeEdges = getRank(msDG, node, True)
                domination_centrality = 2 / (1 / rankControl + 1 / rankObserve)
                msDG = msDG.reverse(copy=True)  # reverse
                # formatting data
                nodeInfo.append((node, "{:.3f}".format(domination_centrality),
                                 rankControl, controlNodes, controlEdges,
                                 rankObserve, observeNodes, observeEdges))
            else:  # node in ms, and is a stem node
                # calculating controllable subspace
                rankControl, controlNodes, controlEdges = getRank(tempCsMsDG, node, False)
                # calculating observable subspace
                rankObserve, observeNodes, observeEdges = getRank(tempOsMsDG, node, True)
                domination_centrality = 2 / (1 / rankControl + 1 / rankObserve)
                # formatting data
                nodeInfo.append((node, "{:.3f}".format(domination_centrality),
                                 rankControl, controlNodes, controlEdges,
                                 rankObserve, observeNodes, observeEdges))

        # store
        # ms count +1
        count += 1

        output_path = os.path.join(tmp_folder, f'./ms/P{i}_ms/P{i}_ms_{count}.txt')
        file = open(output_path, 'w')
        strr = ''
        for data in nodeInfo:
            strr = str(data[0]) + '\t' + str(data[1]) + '\t' + str(data[2]) + '\t' + str(data[3]) + '\t' + str(
                data[4]) + '\t' + str(data[5]) + '\t' + str(data[6]) + '\t' + str(data[7])
            file.write(strr)
            file.write('\n')
        file.close()

print('[7 / 10] Filtering nodes base on control centrality...')

# if not os.path.exists(os.path.join(tmp_folder, 'P1_ms_final.txt')):
core = list(
    pd.read_table(os.path.join(current_directory, f'./data/{args.trait}/omnigenic neighborhood/core.txt'),
                  header=None)[
        0])
per = list(
    pd.read_table(os.path.join(current_directory, f'./data/{args.trait}/omnigenic neighborhood/periphery.txt'),
                  header=None)[0])
per = list(set(per).union(set(core)))

for i in range(1, index + 1):
    ms_path = os.path.join(tmp_folder, f'./ms/P{i}_ms/')
    msFileList = os.listdir(ms_path)

    allMsDC = []
    pbar = tqdm(msFileList, desc=f'Processing P{i}')
    for msFile in pbar:
        msDC = []
        msInfo = open(ms_path + msFile, 'r').readlines()
        for ms in msInfo:
            msDC.append(ms.split('\t')[1])
        allMsDC.append(msDC)
    pbar.close()
    nodeInfoFinal = []
    nodeIndex = 0
    finalIndex = []
    while len(finalIndex) < len(allMsDC[0]):
        tempDC = -1
        numDC = []
        finalFlag = -1
        for msDC in allMsDC:
            numDC.append(msDC[nodeIndex])
        finalIndex.append(numDC.index(max(numDC)))
        nodeIndex += 1

    nodeIndex = 0
    for index in finalIndex:
        msInfo = open(ms_path + msFileList[index], 'r').readlines()
        nodeInfoFinal.append(msInfo[nodeIndex].split('\t'))
        nodeIndex += 1

    output_path = os.path.join(tmp_folder, f'P{i}_ms_final.txt')
    file = open(output_path, 'w')
    strr = ''
    for data in nodeInfoFinal:
        strr = str(data[0]) + '\t' + str(data[1]) + '\t' + str(data[2]) + '\t' + str(data[3]) + '\t' + str(
            data[4]) + '\t' + str(data[5]) + '\t' + str(data[6]) + '\t' + str(data[7])
        file.write(strr)
    #     file.write('\n')
    file.close()

# optimal control structure
print('[8 / 10] Constructing optimal control structure...')
# overall
if not os.path.exists(os.path.join(tmp_folder, './edges/')):
    os.mkdir(os.path.join(tmp_folder, './edges'))

coreNodes = list(pd.read_table(os.path.join(tmp_folder, "P1_coreNodes.txt"), header=None)[0])
C1Nodes = list(pd.read_table(os.path.join(tmp_folder, "P1_neighborNodes.txt"), header=None)[0])
C2Nodes = list(pd.read_table(os.path.join(tmp_folder, "P2_neighborNodes.txt"), header=None)[0])
C3Nodes = list(pd.read_table(os.path.join(tmp_folder, "P3_neighborNodes.txt"), header=None)[0])
C4Nodes = list(pd.read_table(os.path.join(tmp_folder, "P4_neighborNodes.txt"), header=None)[0])

# P1
df = pd.read_table(os.path.join(tmp_folder, 'P1_ms_final.txt'), header=None)
controlEdgesStr = list(df[4])
observeEdgesStr = list(df[7])
nodeInfo = list(df[0])
allEdges = set()
for i in range(len(nodeInfo)):
    controlEdges = readStrToList(controlEdgesStr[i])
    observeEdges = readStrToList(observeEdgesStr[i])
    allEdges = allEdges | set(controlEdges) | set(observeEdges)

COToCO, C1ToCO, noControl, noObserve, others = set(), set(), 0, 0, 0
for i in range(len(nodeInfo)):
    controlEdges = readStrToList(controlEdgesStr[i])
    if len(controlEdges) == 0:
        noControl += 1
    observeEdges = readStrToList(observeEdgesStr[i])
    if len(observeEdges) == 0:
        noObserve += 1
    for edge in set(controlEdges + observeEdges):
        if edge[0] in coreNodes and edge[1] in coreNodes:
            COToCO.add(edge)
        elif edge[0] in C1Nodes and edge[1] in coreNodes:
            C1ToCO.add(edge)
        else:
            others += 1

output_path = os.path.join(tmp_folder, './edges/COToCO_edges.txt')
file = open(output_path, 'w')
strr = ''
for data in COToCO:
    strr = str(data[0]) + '\t' + str('COCO') + '\t' + str(data[1])
    file.write(strr)
    file.write('\n')
file.close()

output_path = os.path.join(tmp_folder, './edges/C1ToCO_edges.txt')
file = open(output_path, 'w')
strr = ''
for data in C1ToCO:
    strr = str(data[0]) + '\t' + str('C1CO') + '\t' + str(data[1])
    file.write(strr)
    file.write('\n')
file.close()

df = pd.read_table(os.path.join(tmp_folder, 'P1_non_binary.txt'), header=None)
source = list(df[0])
target = list(df[1])

output_path = os.path.join(tmp_folder, 'P1_allEdges.txt')
coco_none, coco, c1co_none, c1co, other = 0, 0, 0, 0, 0
file = open(output_path, 'w')
strr = ''
for i in range(len(source)):
    if (source[i], target[i]) in COToCO:
        strr = str(source[i]) + '\t' + str('COCO') + '\t' + str(target[i])
        coco += 1
    elif (source[i], target[i]) in C1ToCO:
        strr = str(source[i]) + '\t' + str('C1CO') + '\t' + str(target[i])
        c1co += 1
    elif source[i] in coreNodes and target[i] in coreNodes:
        strr = str(source[i]) + '\t' + str('COCO-NONE') + '\t' + str(target[i])
        coco_none += 1
    elif source[i] in C1Nodes and target[i] in coreNodes:
        strr = str(source[i]) + '\t' + str('C1CO-NONE') + '\t' + str(target[i])
        c1co_none += 1
    else:
        other += 1
    file.write(strr)
    file.write('\n')
file.close()

# P2
df = pd.read_table(os.path.join(tmp_folder, 'P2_ms_final.txt'), header=None)
controlEdgesStr = list(df[4])
observeEdgesStr = list(df[7])
nodeInfo = list(df[0])
matchEdges = set()
for i in range(len(nodeInfo)):
    controlEdges = readStrToList(controlEdgesStr[i])
    observeEdges = readStrToList(observeEdgesStr[i])
    matchEdges = matchEdges | set(controlEdges) | set(observeEdges)
COToC1, C1ToC1, C2ToC1, noControl, noObserve, others = set(), set(), set(), 0, 0, 0
for i in range(len(nodeInfo)):
    controlEdges = readStrToList(controlEdgesStr[i])
    if len(controlEdges) == 0:
        noControl += 1
    observeEdges = readStrToList(observeEdgesStr[i])
    if len(observeEdges) == 0:
        noObserve += 1
    for edge in set(controlEdges + observeEdges):
        if edge[0] in coreNodes and edge[1] in C1Nodes:
            COToC1.add(edge)
        elif edge[0] in C1Nodes and edge[1] in C1Nodes:
            C1ToC1.add(edge)
        elif edge[0] in C2Nodes and edge[1] in C1Nodes:
            C2ToC1.add(edge)
        else:
            others += 1

output_path = os.path.join(tmp_folder, './edges/COToC1_edges.txt')
file = open(output_path, 'w')
strr = ''
for data in COToC1:
    strr = str(data[0]) + '\t' + str('COC1') + '\t' + str(data[1])
    file.write(strr)
    file.write('\n')
file.close()
output_path = os.path.join(tmp_folder, './edges/C1ToC1_edges.txt')

file = open(output_path, 'w')
strr = ''
for data in C1ToC1:
    strr = str(data[0]) + '\t' + str('C1C1') + '\t' + str(data[1])
    file.write(strr)
    file.write('\n')
file.close()

output_path = os.path.join(tmp_folder, './edges/C2ToC1_edges.txt')
file = open(output_path, 'w')
strr = ''
for data in C2ToC1:
    strr = str(data[0]) + '\t' + str('C2C1') + '\t' + str(data[1])
    file.write(strr)
    file.write('\n')
file.close()

df = pd.read_table(os.path.join(tmp_folder, 'P2_non_binary_selected.txt'), header=None)

source = list(df[0])
target = list(df[1])
output_path = os.path.join(tmp_folder, 'P2_allEdges.txt')
coc1_none, coc1, c1c1_none, c1c1, c2c1_none, c2c1, other = 0, 0, 0, 0, 0, 0, 0
file = open(output_path, 'w')
strr = ''
for i in range(len(source)):
    if (source[i], target[i]) in COToC1:
        strr = str(source[i]) + '\t' + str('COC1') + '\t' + str(target[i])
        coc1 += 1
    elif (source[i], target[i]) in C1ToC1:
        strr = str(source[i]) + '\t' + str('C1C1') + '\t' + str(target[i])
        c1c1 += 1
    elif (source[i], target[i]) in C2ToC1:
        strr = str(source[i]) + '\t' + str('C2C1') + '\t' + str(target[i])
        c2c1 += 1
    elif source[i] in coreNodes and target[i] in C1Nodes:
        strr = str(source[i]) + '\t' + str('COC1-NONE') + '\t' + str(target[i])
        coc1_none += 1
    elif source[i] in C1Nodes and target[i] in C1Nodes:
        strr = str(source[i]) + '\t' + str('C1C1-NONE') + '\t' + str(target[i])
        c1c1_none += 1
    elif source[i] in C2Nodes and target[i] in C1Nodes:
        strr = str(source[i]) + '\t' + str('C2C1-NONE') + '\t' + str(target[i])
        c2c1_none += 1
    else:
        other += 1
    file.write(strr)
    file.write('\n')
file.close()

# join C1 to C2
file = open(output_path, 'a')
c1_path = output_path = os.path.join(tmp_folder, 'P1_allEdges.txt')
file_c1 = open(c1_path, 'r')
for line in file_c1:
    file.write(line)
#     file.write('\n')
file_c1.close()
file.close()

# P3
df = pd.read_table(os.path.join(tmp_folder, 'P3_ms_final.txt'), header=None)

controlEdgesStr = list(df[4])
observeEdgesStr = list(df[7])
nodeInfo = list(df[0])
matchEdges = set()
for i in range(len(nodeInfo)):
    controlEdges = readStrToList(controlEdgesStr[i])
    observeEdges = readStrToList(observeEdgesStr[i])
    matchEdges = matchEdges | set(controlEdges) | set(observeEdges)
COToC2, C1ToC2, C2ToC2, C3ToC2, noControl, noObserve, others = set(), set(), set(), set(), 0, 0, 0
for i in range(len(nodeInfo)):
    controlEdges = readStrToList(controlEdgesStr[i])
    if len(controlEdges) == 0:
        noControl += 1
    observeEdges = readStrToList(observeEdgesStr[i])
    if len(observeEdges) == 0:
        noObserve += 1
    for edge in set(controlEdges + observeEdges):
        if edge[0] in coreNodes and edge[1] in C2Nodes:
            COToC2.add(edge)
        elif edge[0] in C1Nodes and edge[1] in C2Nodes:
            C1ToC2.add(edge)
        elif edge[0] in C2Nodes and edge[1] in C2Nodes:
            C2ToC2.add(edge)
        elif edge[0] in C3Nodes and edge[1] in C2Nodes:
            C3ToC2.add(edge)
        else:
            others += 1
output_path = os.path.join(tmp_folder, './edges/COToC2_edges.txt')

file = open(output_path, 'w')
strr = ''
for data in COToC2:
    strr = str(data[0]) + '\t' + str('COC2') + '\t' + str(data[1])
    file.write(strr)
    file.write('\n')
file.close()
output_path = os.path.join(tmp_folder, './edges/C1ToC2_edges.txt')
file = open(output_path, 'w')
strr = ''
for data in C1ToC2:
    strr = str(data[0]) + '\t' + str('C1C2') + '\t' + str(data[1])
    file.write(strr)
    file.write('\n')
file.close()
output_path = os.path.join(tmp_folder, './edges/C2ToC2_edges.txt')
file = open(output_path, 'w')
strr = ''
for data in C2ToC2:
    strr = str(data[0]) + '\t' + str('C2C2') + '\t' + str(data[1])
    file.write(strr)
    file.write('\n')
file.close()
output_path = os.path.join(tmp_folder, './edges/C3ToC2_edges.txt')
file = open(output_path, 'w')
strr = ''
for data in C3ToC2:
    strr = str(data[0]) + '\t' + str('C3C2') + '\t' + str(data[1])
    file.write(strr)
    file.write('\n')
file.close()
df = pd.read_table(os.path.join(tmp_folder, 'P3_non_binary_selected.txt'), header=None)
source = list(df[0])
target = list(df[1])

output_path = os.path.join(tmp_folder, 'P3_allEdges.txt')
coc2_none, coc2, c1c2_none, c1c2, c2c2_none, c2c2, c3c2_none, c3c2, other = 0, 0, 0, 0, 0, 0, 0, 0, 0
file = open(output_path, 'w')
strr = ''
for i in range(len(source)):
    if (source[i], target[i]) in COToC2:
        strr = str(source[i]) + '\t' + str('COC2') + '\t' + str(target[i])
        coc2 += 1
    elif (source[i], target[i]) in C1ToC2:
        strr = str(source[i]) + '\t' + str('C1C2') + '\t' + str(target[i])
        c1c2 += 1
    elif (source[i], target[i]) in C2ToC2:
        strr = str(source[i]) + '\t' + str('C2C2') + '\t' + str(target[i])
        c2c2 += 1
    elif (source[i], target[i]) in C3ToC2:
        strr = str(source[i]) + '\t' + str('C3C2') + '\t' + str(target[i])
        c3c2 += 1
    elif source[i] in coreNodes and target[i] in C2Nodes:
        strr = str(source[i]) + '\t' + str('COC2-NONE') + '\t' + str(target[i])
        coc2_none += 1
    elif source[i] in C1Nodes and target[i] in C2Nodes:
        strr = str(source[i]) + '\t' + str('C1C2-NONE') + '\t' + str(target[i])
        c1c2_none += 1
    elif source[i] in C2Nodes and target[i] in C2Nodes:
        strr = str(source[i]) + '\t' + str('C2C2-NONE') + '\t' + str(target[i])
        c2c2_none += 1
    elif source[i] in C3Nodes and target[i] in C2Nodes:
        strr = str(source[i]) + '\t' + str('C3C2-NONE') + '\t' + str(target[i])
        c3c2_none += 1
    else:
        other += 1
    file.write(strr)
    file.write('\n')
file.close()

# Join P2 to P3
file = open(output_path, 'a')
# c2_path = "F:/研究生/工作汇报/数据/asthma网络/分层控制结构/画图/asthma_C2/asthma_C2_allEdges.txt"
c2_path = os.path.join(tmp_folder, 'P2_allEdges.txt')
file_c2 = open(c2_path, 'r')
for line in file_c2:
    file.write(line)
#     file.write('\n')
file_c2.close()
file.close()

# P4
df = pd.read_table(os.path.join(tmp_folder, 'P4_ms_final.txt'), header=None)

controlEdgesStr = list(df[4])
observeEdgesStr = list(df[7])
nodeInfo = list(df[0])
matchEdges = set()
for i in range(len(nodeInfo)):
    controlEdges = readStrToList(controlEdgesStr[i])
    observeEdges = readStrToList(observeEdgesStr[i])
    matchEdges = matchEdges | set(controlEdges) | set(observeEdges)
COToC3, C1ToC3, C2ToC3, C3ToC3, C4ToC3, C4ToC4, noControl, noObserve, others = set(), set(), set(), set(), set(), set(), 0, 0, 0
for i in range(len(nodeInfo)):
    controlEdges = readStrToList(controlEdgesStr[i])
    if len(controlEdges) == 0:
        noControl += 1
    observeEdges = readStrToList(observeEdgesStr[i])
    if len(observeEdges) == 0:
        noObserve += 1
    for edge in set(controlEdges + observeEdges):
        if edge[0] in coreNodes and edge[1] in C3Nodes:
            COToC3.add(edge)
        elif edge[0] in C1Nodes and edge[1] in C3Nodes:
            C1ToC3.add(edge)
        elif edge[0] in C2Nodes and edge[1] in C3Nodes:
            C2ToC3.add(edge)
        elif edge[0] in C3Nodes and edge[1] in C3Nodes:
            C3ToC3.add(edge)
        elif edge[0] in C4Nodes and edge[1] in C3Nodes:
            C4ToC3.add(edge)
        elif edge[0] in C4Nodes and edge[1] in C4Nodes:
            C4ToC4.add(edge)
        else:
            others += 1

output_path = os.path.join(tmp_folder, './edges/COToC3_edges.txt')
file = open(output_path, 'w')
strr = ''
for data in COToC3:
    strr = str(data[0]) + '\t' + str('COC3') + '\t' + str(data[1])
    file.write(strr)
    file.write('\n')
file.close()
output_path = os.path.join(tmp_folder, './edges/C1ToC3_edges.txt')
file = open(output_path, 'w')
strr = ''
for data in C1ToC3:
    strr = str(data[0]) + '\t' + str('C1C3') + '\t' + str(data[1])
    file.write(strr)
    file.write('\n')
file.close()
output_path = os.path.join(tmp_folder, './edges/C2ToC3_edges.txt')
file = open(output_path, 'w')
strr = ''
for data in C2ToC3:
    strr = str(data[0]) + '\t' + str('C2C3') + '\t' + str(data[1])
    file.write(strr)
    file.write('\n')
file.close()
output_path = os.path.join(tmp_folder, './edges/C3ToC3_edges.txt')
file = open(output_path, 'w')
strr = ''
for data in C3ToC3:
    strr = str(data[0]) + '\t' + str('C3C3') + '\t' + str(data[1])
    file.write(strr)
    file.write('\n')
file.close()
output_path = os.path.join(tmp_folder, './edges/C4ToC3_edges.txt')
file = open(output_path, 'w')
strr = ''
for data in C4ToC3:
    strr = str(data[0]) + '\t' + str('C4C3') + '\t' + str(data[1])
    file.write(strr)
    file.write('\n')
file.close()
output_path = os.path.join(tmp_folder, './edges/C4ToC4_edges.txt')
file = open(output_path, 'w')
strr = ''
for data in C4ToC4:
    strr = str(data[0]) + '\t' + str('C4C4') + '\t' + str(data[1])
    file.write(strr)
    file.write('\n')
file.close()

df = pd.read_table(os.path.join(tmp_folder, 'P4_non_binary_selected.txt'), header=None)

source = list(df[0])
target = list(df[1])

output_path = os.path.join(tmp_folder, 'P4_allEdges.txt')
coc3_none, coc3, c1c3_none, c1c3, c2c3_none, c2c3, c3c3_none, c3c3, c4c3_none, c4c3, c4c4_none, c4c4, other = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
file = open(output_path, 'w')
strr = ''
for i in range(len(source)):
    if (source[i], target[i]) in COToC3:
        strr = str(source[i]) + '\t' + str('COC3') + '\t' + str(target[i])
        coc3 += 1
    elif (source[i], target[i]) in C1ToC3:
        strr = str(source[i]) + '\t' + str('C1C3') + '\t' + str(target[i])
        c1c3 += 1
    elif (source[i], target[i]) in C2ToC3:
        strr = str(source[i]) + '\t' + str('C2C3') + '\t' + str(target[i])
        c2c3 += 1
    elif (source[i], target[i]) in C3ToC3:
        strr = str(source[i]) + '\t' + str('C3C3') + '\t' + str(target[i])
        c3c3 += 1
    elif (source[i], target[i]) in C4ToC3:
        strr = str(source[i]) + '\t' + str('C4C3') + '\t' + str(target[i])
        c4c3 += 1
    elif (source[i], target[i]) in C4ToC4:
        strr = str(source[i]) + '\t' + str('C4C4') + '\t' + str(target[i])
        c4c4 += 1
    elif source[i] in coreNodes and target[i] in C3Nodes:
        strr = str(source[i]) + '\t' + str('COC3-NONE') + '\t' + str(target[i])
        coc3_none += 1
    elif source[i] in C1Nodes and target[i] in C3Nodes:
        strr = str(source[i]) + '\t' + str('C1C3-NONE') + '\t' + str(target[i])
        c1c3_none += 1
    elif source[i] in C2Nodes and target[i] in C3Nodes:
        strr = str(source[i]) + '\t' + str('C2C3-NONE') + '\t' + str(target[i])
        c2c3_none += 1
    elif source[i] in C3Nodes and target[i] in C3Nodes:
        strr = str(source[i]) + '\t' + str('C3C3-NONE') + '\t' + str(target[i])
        c3c3_none += 1
    elif source[i] in C4Nodes and target[i] in C3Nodes:
        strr = str(source[i]) + '\t' + str('C4C3-NONE') + '\t' + str(target[i])
        c4c3_none += 1
    elif source[i] in C4Nodes and target[i] in C4Nodes:
        strr = str(source[i]) + '\t' + str('C4C4-NONE') + '\t' + str(target[i])
        c4c4_none += 1
    else:
        other += 1
    file.write(strr)
    file.write('\n')
file.close()

file = open(output_path, 'a')
c3_path = os.path.join(tmp_folder, 'P3_allEdges.txt')
file_c3 = open(c3_path, 'r')
for line in file_c3:
    file.write(line)
#     file.write('\n')
file_c3.close()
file.close()

# Edge criticality
print('[9 / 10] Calculating edge criticality...')

coreNodes = list(pd.read_table(os.path.join(tmp_folder, "P1_coreNodes.txt"), header=None)[0])
P1Nodes = list(pd.read_table(os.path.join(tmp_folder, "P1_neighborNodes.txt"), header=None)[0])
P2Nodes = list(pd.read_table(os.path.join(tmp_folder, "P2_neighborNodes.txt"), header=None)[0])
P3Nodes = list(pd.read_table(os.path.join(tmp_folder, "P3_neighborNodes.txt"), header=None)[0])
P4Nodes = list(pd.read_table(os.path.join(tmp_folder, "P4_neighborNodes.txt"), header=None)[0])

allNodes = coreNodes + P1Nodes + P2Nodes + P3Nodes + P4Nodes

nodeProbability = {}

pro_path = os.path.join(current_directory, f'./data/{args.trait}/genetic perturbation values.txt')
pro_file = open(pro_path, 'r')
for line in pro_file:
    pas = line.strip().split()
    if pas[0] in allNodes:
        nodeProbability.update({pas[0]: float(pas[1])})
pro_file.close()

originNetwork = nx.DiGraph()
path = os.path.join(tmp_folder, 'before_binary.txt')
file = open(path, 'r')
for line in file:
    pas = line.strip().split()
    originNetwork.add_edge(pas[0], pas[1])
file.close()

controlNetwork = nx.DiGraph()
controlEdges = []
C4_path = os.path.join(tmp_folder, 'P4_allEdges.txt')
C4_file = open(C4_path, 'r')
for line in C4_file:
    pas = line.strip().split()
    if 'NONE' in pas[1]:
        continue
    controlNetwork.add_edge(pas[0], pas[2])
    controlEdges.append((pas[0], pas[1], pas[2]))

output_path = os.path.join(tmp_folder, 'control_network.txt')
file = open(output_path, 'w')
strr = ''
for data in controlEdges:
    strr = str(data[0]) + '\t' + str(data[1]) + '\t' + str(data[2])
    file.write(strr)
    file.write('\n')
file.close()

nodeLayer = {}
for node in coreNodes:
    nodeLayer.update({node: 0})
for node in P1Nodes:
    nodeLayer.update({node: 1})
for node in P2Nodes:
    nodeLayer.update({node: 2})
for node in P3Nodes:
    nodeLayer.update({node: 3})
for node in P4Nodes:
    nodeLayer.update({node: 4})
layerList = [coreNodes, P1Nodes, P2Nodes, P3Nodes, P4Nodes]
deleteNetwork = nx.DiGraph()
deleteNetwork.add_edges_from(list(controlNetwork.edges()))
controlEdges = list(controlNetwork.edges())
index = 1
edgeScore = {}
pbar = tqdm(controlEdges, desc='Processing')
for edge in pbar:
    index += 1

    observeNum = len(nx.ancestors(deleteNetwork, edge[0]))

    controlBeforeDelete = getNodeControlScore(deleteNetwork, edge[0], nodeLayer, layerList, nodeProbability)

    deleteNetwork.remove_edge(edge[0], edge[1])

    controlAfterDelete = getNodeControlScore(deleteNetwork, edge[0], nodeLayer, layerList, nodeProbability)
    edgeControlScore = observeNum * abs(float(controlAfterDelete) - float(controlBeforeDelete))

    deleteNetwork.add_edge(edge[0], edge[1])

    edgeObserveScore = 0
    edgeScore.update({edge: (edgeControlScore, edgeObserveScore)})
pbar.close()
if not os.path.exists(os.path.join(current_directory, f'./output/{args.trait}/')):
    os.mkdir(os.path.join(current_directory, f'./output/{args.trait}/'))

output_path = os.path.join(current_directory, f'./output/{args.trait}/edge_score.txt')
file = open(output_path, 'w')
strr = ''
for key, value in edgeScore.items():
    strr = str(key[0]) + '\t' + str(key[1]) + '\t' + str(value[0]) + '\t' + str(value[1])
    file.write(strr)
    file.write('\n')
file.close()

# edge removal sequence
print('[10 / 10] Calculating descend of control ability...')

ascend_score = sorted(edgeScore.items(), key=lambda x: x[1][0], reverse=False)

deleteNetwork = nx.DiGraph()

deleteNetwork.add_edges_from(list(controlNetwork.edges()))

allNodes = set(controlNetwork.nodes())
ascend_change = []
deleteNodes = set()
i = 1
C4C3, C3C2, C2C1, C1CO = 0, 0, 0, 0
for data in ascend_score:
    i += 1
    if data[0][0] in P4Nodes and data[0][1] in P3Nodes:
        C4C3 += 1
    if data[0][0] in P3Nodes and data[0][1] in P2Nodes:
        C3C2 += 1
    if data[0][0] in P2Nodes and data[0][1] in P1Nodes:
        C2C1 += 1
    if data[0][0] in P1Nodes and data[0][1] in coreNodes:
        C1CO += 1
    deleteNetwork.remove_edge(data[0][0], data[0][1])
    num = getControlNodeNum(set(allNodes), set(coreNodes), deleteNodes, deleteNetwork)  # 查看可控核心区域的平均节点数
    ascend_change.append(
        (num, len(deleteNetwork.nodes()), str(C4C3), str(C3C2), str(C2C1), str(C1CO), len(deleteNodes)))

output_path = os.path.join(current_directory, f'./output/{args.trait}/delete_control_score.txt')

file = open(output_path, 'w')
for score in ascend_change:
    strr = str(score[0]) + '\t' + str(score[1]) + '\t' + str(score[2]) + '\t' + str(score[3]) + '\t' + str(
        score[4]) + '\t' + str(score[5]) + '\t' + str(score[6])
    file.write(strr)
    file.write('\n')
file.close()

print("Done.")
