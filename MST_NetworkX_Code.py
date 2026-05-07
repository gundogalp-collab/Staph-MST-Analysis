import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from itertools import combinations
import numpy as np
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import re
import matplotlib.patheffects as PathEffects

# =========================
# 1️⃣ File Reading
# =========================
df = pd.read_excel("Supplementary File.xlsx", header=2)
df.columns = df.columns.str.strip()

spa_col = "spa-type"
repeat_col = "Spa Repeats"
source_col = "Sample Source"

df = df[[spa_col, repeat_col, source_col]].dropna()

# =========================
# 2️⃣ Filter & Feedback
# =========================

# ==========================================
# 2. Filtering 
# ==========================================
def get_clean_tokens(text):
    if not isinstance(text, str): return []
    clean = re.sub(r'[-,]', ' ', text)
    return [t for t in clean.split() if t]

df['repeat_tokens'] = df[repeat_col].apply(get_clean_tokens)
df['repeat_len'] = df['repeat_tokens'].apply(len)

df_filtered = df[df[spa_col].astype(str) != "NT"].copy()

df = df_filtered

# =========================
# 3️⃣ Group
# =========================
grouped = df.groupby(spa_col)
spa_types = []
rep_dict = {} 
counts = {}

for name, group in grouped:
    spa_types.append(name)
    rep_dict[name] = group['repeat_tokens'].iloc[0]
    counts[name] = len(group) 

token_dict = rep_dict

# =========================
# 4️⃣ Edit Distance (Levenshtein) 
# =========================
def token_levenshtein(a, b):
    m, n = len(a), len(b)
    dp = [[0]*(n+1) for _ in range(m+1)]
    for i in range(m+1): dp[i][0] = i
    for j in range(n+1): dp[0][j] = j
    for i in range(1, m+1):
        for j in range(1, n+1):
            cost = 0 if a[i-1] == b[j-1] else 1
            dp[i][j] = min(dp[i-1][j]+1, dp[i][j-1]+1, dp[i-1][j-1]+cost)
    return dp[m][n]

G = nx.Graph()
for spa in spa_types: G.add_node(spa)

for a, b in combinations(spa_types, 2):
    d = token_levenshtein(token_dict[a], token_dict[b])
    G.add_edge(a, b, weight=d)

if len(G.nodes) > 0:
    # Minimum Spanning Tree
    mst = nx.minimum_spanning_tree(G, weight="weight")
    
    # =========================
    # 5️⃣ Visualizing
    # =========================
    source_distribution = df.groupby([spa_col, source_col]).size().unstack(fill_value=0).to_dict("index")
    sources = df[source_col].unique()
    base_colors = plt.cm.tab10.colors


    color_map = {
        "Hospital": "#D90429",         
        "Dairy Farm": "#0077B6",       
        "Dairy Plant": "#00B4D8",      
        "Livestock farm": "#2A9D8F",   
        "Slaughterhouse": "#E76F51",   
        "Meat Plant": "#6A4C93"        
        }

    node_radii = {n: 0.03 + (np.sqrt(counts[n]) * 0.015) for n in mst.nodes()}
    
    pos = nx.spring_layout(mst, k=3.5, iterations=1000, seed=42)

    def prevent_overlap(pos, radii, padding=0.06, iterations=50):
        pos_new = pos.copy()
        keys = list(pos_new.keys())
        for _ in range(iterations):
            moved = False
            for i in range(len(keys)):
                for j in range(i + 1, len(keys)):
                    n1, n2 = keys[i], keys[j]
                    x1, y1 = pos_new[n1]
                    x2, y2 = pos_new[n2]
                    dx, dy = x1 - x2, y1 - y2
                    dist = np.sqrt(dx*dx + dy*dy)
                    min_dist = radii[n1] + radii[n2] + padding
                    if dist < min_dist:
                        moved = True
                        if dist == 0: dx, dy = 0.01, 0; dist = 0.01
                        overlap = min_dist - dist
                        push = overlap / 2
                        pos_new[n1] = (x1 + (dx/dist)*push, y1 + (dy/dist)*push)
                        pos_new[n2] = (x2 - (dx/dist)*push, y2 - (dy/dist)*push)
            if not moved: break
        return pos_new

    pos = prevent_overlap(pos, node_radii, padding=0.02)

    # =========================
    # 6️⃣ Draw
    # =========================
    plt.figure(figsize=(18, 16))
    plt.gca().set_aspect("equal")

    edges_thick =  [(u,v) for u,v,d in mst.edges(data=True) if d['weight'] <= 2.2]
    edges_thinner =[(u,v) for u,v,d in mst.edges(data=True) if 2.2 < d['weight'] <= 3.4]
    edges_thinnest=[(u,v) for u,v,d in mst.edges(data=True) if 3.4 < d['weight'] <= 4.59]
    edges_dashed = [(u,v) for u,v,d in mst.edges(data=True) if 4.59 < d['weight'] <= 5.8]
    edges_dotted = [(u,v) for u,v,d in mst.edges(data=True) if d['weight'] > 5.8]

    nx.draw_networkx_edges(mst, pos, edgelist=edges_thick, width=3.5, edge_color="black", style="solid", alpha=0.8)
    nx.draw_networkx_edges(mst, pos, edgelist=edges_thinner, width=2.0, edge_color="black", style="solid", alpha=0.8)
    nx.draw_networkx_edges(mst, pos, edgelist=edges_thinnest, width=1.0, edge_color="black", style="solid", alpha=0.8)
    nx.draw_networkx_edges(mst, pos, edgelist=edges_dashed, width=1.2, edge_color="black", style="dashed", alpha=0.8)
    nx.draw_networkx_edges(mst, pos, edgelist=edges_dotted, width=1.2, edge_color="black", style="dotted", alpha=0.8)

    for node in mst.nodes():
        x, y = pos[node]
        sizes = list(source_distribution[node].values())
        labels = list(source_distribution[node].keys())
        total = sum(sizes)
        radius = node_radii[node]
        
        start = 0
        for frac, label in zip([s/total for s in sizes], labels):
            theta = frac * 360
            wedge = mpatches.Wedge((x, y), radius, start, start+theta, facecolor=color_map[label], edgecolor="black", linewidth=0.5, zorder=10)
            plt.gca().add_patch(wedge)
            start += theta

    for node in mst.nodes():
        x, y = pos[node]
        font_size = 6 if counts[node] < 5 else 8
        
        # Siyah yazı, ama etrafında beyaz parlama efekti var
        txt = plt.text(x, y, node, ha='center', va='center', fontsize=font_size, 
                       fontweight='bold', zorder=10, color='black')
        txt.set_path_effects([PathEffects.withStroke(linewidth=2.5, foreground='w')])
    # =========================
    # 7️⃣ Legend
    # =========================

    source_patches = [mpatches.Patch(color=color_map[s], label=s) for s in sources]
    first_legend = plt.legend(handles=source_patches, title="Sample Source", loc="upper left", bbox_to_anchor=(-0.25,1), frameon=True)
    plt.gca().add_artist(first_legend)

    line_lines = [
        Line2D([0], [0], color='black', lw=3.5, label='Distance ≤ 2.2'),
        Line2D([0], [0], color='black', lw=2.0, label='2.2 < Distance ≤ 3.4'),
        Line2D([0], [0], color='black', lw=1.0, label='3.4 < Distance ≤ 4.59'),
        Line2D([0], [0], color='black', ls='dashed', label='4.59 < Distance ≤ 5.8'),
        Line2D([0], [0], color='black', ls='dotted', label='Distance > 5.8')
    ]
    plt.legend(handles=line_lines, title="Genetic Distance", loc="lower right", frameon=True)

    #plt.title(f"Minimum Spanning Tree of spa Types\n(Node size ∝ Frequency; Filter: Repeats ≥ 5, No NT)", fontsize=16, fontweight='bold')
    plt.axis("off")
    plt.tight_layout()
    plt.savefig("spa_MST_Strict_Sized.png", dpi=600, bbox_inches="tight")
    plt.show()

else:
    print("No data found")
    print("Check repeat datas")
