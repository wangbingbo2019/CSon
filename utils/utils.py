import random
from queue import Queue
import networkx as nx
from matplotlib.colors import Normalize


def getNeighborsFromNodes(network, nodes):
    network = network.reverse(copy=True)  # reverse
    neighbors = set()
    neighborEdges = []
    for node in nodes:
        neighbors = neighbors | set(network.neighbors(node))
        for neighbor in set(network.neighbors(node)):
            neighborEdges.append((neighbor, node))
    network = network.reverse(copy=True)  # reverse
    return neighbors, neighborEdges


def MsToList(path):  # convert ms.txt to list
    tempList1 = []
    index = 0
    with open(path, 'r') as f:
        while True:
            line = f.readline()
            if not line:
                break

            index += 1
            if len(tempList1) < 1000 and random.random() >= 0.5:
                #             if len(tempList1) < 7000:
                tempList1.append(line.lstrip('"[').rstrip(']\n"'))
                if len(tempList1) == 1000:
                    break

    tempList2 = []
    for t in range(len(tempList1)):
        tempList2.append([[subitem.split(", ")[0].strip("'"), subitem.split(", ")[1].strip("'")] for subitem in
                          [item.lstrip(', (').lstrip('"(') for item in tempList1[t].split(')')][:-1]])

    msList = []
    for t in range(len(tempList2)):
        msList.append([(item[0][:-4], item[1][:-3]) for item in tempList2[t]])
    return msList


def getRank(G, node, matching, reverse):
    if node is None:
        return
    queue = Queue()
    nodeSet = set()
    rank = 1
    edgeSet = set()
    result = []
    queue.put(node)
    nodeSet.add(node)
    while not queue.empty():
        cur = queue.get()
        for neighbor in G.neighbors(cur):
            if not reverse:
                if ((cur, neighbor) in matching) and (neighbor not in nodeSet):
                    rank = rank + 1
                    nodeSet.add(neighbor)
                    edgeSet.add((cur, neighbor))
                    queue.put(neighbor)
            else:
                if ((neighbor, cur) in matching) and (neighbor not in nodeSet):
                    rank = rank + 1
                    nodeSet.add(neighbor)
                    edgeSet.add((cur, neighbor))
                    queue.put(neighbor)
    result.append(rank)
    result.append(nodeSet)
    result.append(edgeSet)
    return result


def getRank(msDG, node, isReverse):
    if node is None:
        return

    subGraph = nx.subgraph(msDG, list(nx.descendants(msDG, node)) + list([node]))
    if isReverse:
        subGraph = subGraph.reverse(copy=True)

    result = []
    rank = len(subGraph.nodes)
    nodeList = list(subGraph.nodes)
    edgeList = list(subGraph.edges)
    return [rank, nodeList, edgeList]


def readStrToList(edges):
    return [(subitem.split(", ")[0].strip("'"), subitem.split(", ")[1].strip("'")) for subitem in [item.lstrip(', (').lstrip('"(') for item in edges.lstrip('"[').rstrip(']\n"').split(')')][:-1]]


def getNodeControlScore(network, node, nodeLayer, layerList, posteriorInfor):
    sourceLayer = nodeLayer[node]
    controlNodes = set(nx.descendants(network, node))
    controlCore = controlNodes & set(layerList[0])
    controlScore = 0
    for core in controlCore:
        paths = list(nx.all_shortest_paths(network, source=node, target=core))
        number = len(paths)
        length = len(paths[0])
        controlScore += (posteriorInfor[core] * (sourceLayer + 1)/length * 1/number)
    return controlScore


def getNodeObeserveScore(network, node, nodeLayer, layerList, posteriorInfor):
    targetLayer = nodeLayer[node]
    observeNodes = nx.ancestors(network, node)
    observeScore = 0
    for observe in observeNodes:
        sourceLayer = nodeLayer[observe]
        # 获取最短路径
        paths = list(nx.all_shortest_paths(network, source=observe, target=node))
        number = len(paths)
        length = len(paths[0])
        if sourceLayer >= targetLayer:
            observeScore += (posteriorInfor[observe] * (sourceLayer - targetLayer + 1)/length * 1/number)
        else:
            observeScore += (posteriorInfor[observe] * 2/length * 1/number)
    return observeScore


def getControlNodeNum(allNodes, coreNodes, deleteNodes, network):
    num = 0 
    for node in (allNodes - deleteNodes):
        descendantNodes = nx.descendants(network, node) | {node}
        controlCoreNum = len(descendantNodes & coreNodes)
        if controlCoreNum != 0:
            num  += controlCoreNum
        else :
            deleteNodes.add(node)
    return num/len(network.nodes())


def find_max_drop_index(nums):
    if len(nums) < 2:
        return None

    max_diff = nums[0] - nums[1]
    max_index = 0

    for i in range(1, len(nums) - 1):
        current_diff = nums[i] - nums[i + 1]
        if current_diff > max_diff:
            max_diff = current_diff
            max_index = i

    return max_index


class RearKNormalize(Normalize):
    def __init__(self, vmin, vmax, rear_k):
        super().__init__(vmin=vmin, vmax=vmax)
        self.rear_k = rear_k

    def __call__(self, val):
        result = super().__call__(value=val)
        result = result * self.rear_k + (1 - self.rear_k)
        return result

