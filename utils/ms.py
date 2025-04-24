import random
import os
import math
import pandas as pd
from copy import deepcopy
import time
from tqdm import tqdm


def genInitialMaxMatch(edge, nx, ny, mx, my):
    # 使用dfs版本的匈牙利算法产生初始的最大匹配
    visited = {}
    MMS = []
    random.shuffle(nx)
    random.shuffle(ny)  # 随机打乱节点的访问顺序，使得产生的最大匹配是随机的
    for i in nx:
        if mx[i] == -1:  # 对还没有匹配节点的左边的节点
            for key in ny:  # 遍历右边的所有节点
                visited[key] = 0  # 初始时将右边的节点都设置为未访问
            path(edge, mx, my, visited, MMS, i)
    return MMS


def genInitialMaxMatchAdjust(edge, nx, ny, mx, my, MMS, target):
    # 使用dfs版本的匈牙利算法产生初始的最大匹配
    visited = {}
    #     MMS=[]
    random.shuffle(nx)
    random.shuffle(ny)  # 随机打乱节点的访问顺序，使得产生的最大匹配是随机的
    for i in nx:
        if mx[i] == -1:  # 对还没有匹配节点的左边的节点
            for key in ny:  # 遍历右边的所有节点
                visited[key] = 0  # 初始时将右边的节点都设置为未访问
            visited[target] = 1
            path(edge, mx, my, visited, MMS, i)
    return MMS


def path(edge, mx, my, visited, MMS, u):
    # 使用dfs寻找从节点u开始的增广路径，若存在，则返回true，否则返回false
    for v in edge[u]:
        if not visited[v]:  # 存在左边节点到右边节点的边(u,v)且v还没有被从u开始的路径访问
            visited[v] = 1
            if my[v] == -1:  # 如果节点v还没有匹配节点，就直接把u,v设为匹配节点
                mx[u] = v
                my[v] = u
                MMS.append((u, v))
                global linkConnectedToNewMatchedNode
                linkConnectedToNewMatchedNode = (u, v)  # 将in set中新增的匹配节点保存起来，便于之后再根据新增的匹配节点产生新的MMS
                return True
            else:  # 如果v已经有匹配节点u1了
                if path(edge, mx, my, visited, MMS, my[v]):  # 如果重新找到了节点u1的匹配节点，就把u,v设为匹配节点
                    MMS.remove((my[v], v))  # 把v之前的匹配边(u1,v)删除
                    mx[u] = v
                    my[v] = u
                    MMS.append((u, v))
                    return True
    return False


def genMaxMatchByReplacement(edge, ny, mx, my, MMS, j):
    # 在已有MMS的基础上，替换掉j原来的匹配节点，从而得到新的MMS，如果找到了新的MMS，返回in set中新增节点对应的匹配边(已经保存在了全局变量linkConnectedToNewMatchedNode中),否则返回false；新的MMS就是调用方法后的MMS
    visited = {}
    for key in ny:
        visited[key] = 0  # 表示ny中的节点都还没有出现在j开始的增广路径上
    if path(edge, mx, my, visited, MMS, j):
        global linkConnectedToNewMatchedNode
        return linkConnectedToNewMatchedNode
    else:
        return False


def randomly_generate_next_M(edge, ny, mx, my, MMS, i):
    '''
    替换掉当前MMS中in set中的一个节点i，随机产生下一个最大匹配并返回
    '''
    # 删除与节点i相关的所有link
    removeAllLinksConnectedToI(edge, i)
    # 找出与节点i匹配的节点j
    j = my[i]
    # 重置mx,my
    mx[j] = -1
    my[i] = -1
    # 更新MMS
    # print MMS
    MMS.remove((j, i))
    # print(j,i)
    random.shuffle(edge[j])  # 将节点j在in set中的后继节点的顺序打乱，从而随机产生一个新的最大匹配
    # 如果返回False，则没有产生新的最大匹配；否则新的最大匹配就是MMS
    result = genMaxMatchByReplacement(edge, ny, mx, my, MMS, j)  # 新增in set中的匹配节点对应的匹配边
    # 如果返回False，说明无法找到节点j的新的匹配节点，从而产生新的最大匹配，所以返回前须将原最大匹配还原
    if result == False:
        MMS.append((j, i))
        mx[j] = i
        my[i] = j
    return result


def removeAllLinksConnectedToI(edge, i):
    # 删除edge中与节点i相关的所有links
    it = edge.items()
    for (x, ys) in it:
        if i in ys:
            ys.remove(i)


def read_diGraph_edges_from_file(filename):
    '''从文件中读取有向图的所有边，并返回邻接表表示'''
    edge = {}
    f = open(filename, 'r')
    for line in f:
        e = line.strip().split('\t')
        if e[0] == e[1]:
            continue
        if e[0] not in edge:
            edge[e[0]] = []
        edge[e[0]].append(e[1])
    return edge


def get_source_sink_nodes(edge):
    '''根据有向图的邻接表表示，返回所有的源点和汇点'''
    nx = set()
    ny = set()
    for x, ys in edge.items():
        nx.add(x)
        ny = ny | set(ys)
    return list(nx), list(ny)


def enum(digraph_edges_filename, index):
    """ enumerate effective number of max matching from directed graph"""
    edge = read_diGraph_edges_from_file(digraph_edges_filename)
    nx, ny = get_source_sink_nodes(edge)

    nodes_count = len(set(nx) | set(ny))

    # print len(edge),len(nx),len(ny),len(mx),len(my)
    # print len(genInitialMaxMatch(edge,nx,ny,mx,my))
    Ms = set()
    Ms_after = set()
    # consider all edges
    df = pd.read_table(digraph_edges_filename, header=None)
    sourceNode = list(df[0])
    targetNode = list(df[1])

    pbar = tqdm(range(1000), desc=f'Enumerating P{index} Ms', mininterval=0)
    
    for _ in pbar:
        MMS = []
        mx = dict(zip(nx, [-1] * len(nx)))
        my = dict(zip(ny, [-1] * len(ny)))
        choice = random.randint(0, len(sourceNode))
        mx[sourceNode[choice]] = targetNode[choice]
        my[targetNode[choice]] = sourceNode[choice]
        MMS.append((sourceNode[choice], targetNode[choice]))
        M = genInitialMaxMatchAdjust(edge, nx, ny, mx, my, MMS, targetNode[choice])
        Ms.add(str(sorted(M)))

        # randomly choosing node from 'in set' and replace it
        removal = random.choice(M)[1]
        # effective enumerate times
        iterations = int(nodes_count * math.log(nodes_count))

        # outer while for ms enumerate
        for i in range(iterations):
            edge_copy = deepcopy(edge)
            result = randomly_generate_next_M(edge_copy, ny, mx, my, M, removal)

            while result:  # found new ms,
                Ms_after.add(str(sorted(M)))
                removal = result[1]  # 将新加入的被匹配节点替换掉，产生新的最大匹配
                result = randomly_generate_next_M(edge_copy, ny, mx, my, M, removal)  # mx, my correspond to M

            removal = random.choice(M)[1]  # 从当前最大匹配中随机挑选出新的removal节点
    pbar.close()
    Ms = Ms | Ms_after
    return Ms


def save_Ms(Ms, fileout_name):
    f = open(fileout_name, 'w')
    for M in Ms:
        f.write(M.strip('""') + '\n')
    f.close()


def enum_and_save_Ms(tmp_folder, index):
    Ms = enum(os.path.join(tmp_folder, 'P1_binary.txt'), 1)
    save_Ms(Ms, os.path.join(tmp_folder, 'P1_maxmatching.txt'))
    for i in range(2, index + 1):
        Ms = enum(os.path.join(tmp_folder, f'P{i}_binary_selected.txt'), i)
        save_Ms(Ms, os.path.join(tmp_folder, f'P{i}_maxmatching.txt'))
