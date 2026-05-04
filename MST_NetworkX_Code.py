import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from itertools import combinations
import numpy as np
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import re

# =========================
# 1.File Reading
# =========================
df = pd.read_excel("Supplementary File.xlsx", header=2)
df.columns = df.columns.str.strip()

spa_col = "spa-type"
repeat_col = "Spa Repeats"
source_col = "Sample Source"

df = df[[spa_col, repeat_col, source_col]].dropna()

# =========================
# 2.Filtering and report
# =========================

def get_clean_tokens(text):
    if not isinstance(text, str): return []
    clean = re.sub(r'[-,]', ' ', text)
    return [t for t in clean.split() if t]

df['repeat_tokens'] = df[repeat_col].apply(get_clean_tokens)
df['repeat_len'] = df['repeat_tokens'].apply(len)

print("\n" + "="*50)
print("EXCLUSION REPORT")
print("="*50)

# 1. Find NT's
nt_data = df[df[spa_col].astype(str) == "NT"]
if not nt_data.empty:
    print(f"\n1️⃣  [NT] Non-Typeable ones:")
    print(f"    A total of {len(nt_data)} NT's were deleted")
else:
    print("\n1️⃣  NT (Non-Typeable)")

# 2. Find short repeats (<5) 
short_repeats_data = df[
    (df[spa_col].astype(str) != "NT") & 
    (df['repeat_len'] < 5)
]

if not short_repeats_data.empty:
    print(f"\n2️⃣  [SHORT] Short Repeats")
    short_groups = short_repeats_data.groupby(spa_col)
    
    for spa_name, group in short_groups:
        length = group['repeat_len'].iloc[0]
        count = len(group)
        sequence = group[repeat_col].iloc[0]
        print(f"   ❌ {spa_name:<8} | Length: {length} | Isolate count: {count} | Sequence: {sequence}")
else:
    print("\n2️⃣  non short repeats were found")

print("="*50 + "\n")

df_filtered = df[
    (df[spa_col].astype(str) != "NT") & 
    (df['repeat_len'] >= 5) 
].copy()

# General Report
print("--- FINAL STATISTICS ---")
print(f"Total start data : {len(df)}")
print(f"Left after filtering    : {len(df_filtered)}")
print(f"Spa types that entered analyses  : {df_filtered[spa_col].nunique()}")
print("-" * 30)

df = df_filtered

# =========================
# 3️.Groupping
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
# 4️.Edit Distance (Levenshtein) 
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
    # 5. Visualition
    # =========================
    source_distribution = df.groupby([spa_col, source_col]).size().unstack(fill_value=0).to_dict("index")
    sources = df[source_col].unique()
    base_colors = plt.cm.tab10.colors
    color_map = {source: base_colors[i % len(base_colors)] for i, source in enumerate(sources)}
    if "Hospital" in color_map: color_map["Hospital"] = "#FFD700"

    # --- NODE SPACING ---
    node_radii = {n: 0.03 + (np.sqrt(counts[n]) * 0.015) for n in mst.nodes()}
    
    # Layout Calculating
    pos = nx.spring_layout(mst, k=1.5, iterations=1000, seed=42)

    # Overlap preventing
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
    # 6.Drawing
    # =========================
    plt.figure(figsize=(18, 16))

    # Edges
    edges_thick =  [(u,v) for u,v,d in mst.edges(data=True) if d['weight'] <= 2.2]
    edges_thinner =[(u,v) for u,v,d in mst.edges(data=True) if 2.2 < d['weight'] <= 3.4]
    edges_thinnest=[(u,v) for u,v,d in mst.edges(data=True) if 3.4 < d['weight'] <= 4.59]
    edges_dashed = [(u,v) for u,v,d in mst.edges(data=True) if 4.59 < d['weight'] <= 5.8]
    edges_dotted = [(u,v) for u,v,d in mst.edges(data=True) if d['weight'] > 5.8]

    # Line styles
    nx.draw_networkx_edges(mst, pos, edgelist=edges_thick, width=3.5, edge_color="black", style="solid", alpha=0.8)
    nx.draw_networkx_edges(mst, pos, edgelist=edges_thinner, width=2.0, edge_color="black", style="solid", alpha=0.8)
    nx.draw_networkx_edges(mst, pos, edgelist=edges_thinnest, width=1.0, edge_color="black", style="solid", alpha=0.8)
    nx.draw_networkx_edges(mst, pos, edgelist=edges_dashed, width=1.2, edge_color="black", style="dashed", alpha=0.8)
    nx.draw_networkx_edges(mst, pos, edgelist=edges_dotted, width=1.2, edge_color="black", style="dotted", alpha=0.8)

    # Nodes
    for node in mst.nodes():
        x, y = pos[node]
        sizes = list(source_distribution[node].values())
        labels = list(source_distribution[node].keys())
        total = sum(sizes)
        radius = node_radii[node]
        
        start = 0
        for frac, label in zip([s/total for s in sizes], labels):
            theta = frac * 360
            wedge = mpatches.Wedge((x, y), radius, start, start+theta, facecolor=color_map[label], edgecolor="black", linewidth=0.5)
            plt.gca().add_patch(wedge)
            start += theta

    # Etiketler
    for node in mst.nodes():
        x, y = pos[node]
        font_size = 7 if counts[node] < 5 else 9
        plt.text(x, y, node, ha='center', va='center', fontsize=font_size, fontweight='bold', zorder=10, 
                 bbox=dict(facecolor='white', edgecolor='none', alpha=0.6, boxstyle='round,pad=0.1'))

    # =========================
    # 7️.Legend
    # =========================
    
    # Sources
    source_patches = [mpatches.Patch(color=color_map[s], label=s) for s in sources]
    first_legend = plt.legend(handles=source_patches, title="Sample Source", loc="upper left", frameon=True)
    plt.gca().add_artist(first_legend)
    
    # Distances
    line_lines = [
        Line2D([0], [0], color='black', lw=3.5, label='Diff ≤ 2.2'),
        Line2D([0], [0], color='black', lw=2.0, label='2.2 < Diff ≤ 3.4'),
        Line2D([0], [0], color='black', lw=1.0, label='3.4 < Diff ≤ 4.59'),
        Line2D([0], [0], color='black', ls='dashed', label='4.59 < Diff ≤ 5.8'),
        Line2D([0], [0], color='black', ls='dotted', label='Diff > 5.8')
    ]
    plt.legend(handles=line_lines, title="Genetic Distance", loc="lower right", frameon=True)

    plt.title(f"Minimum Spanning Tree of spa Types\n(Node size ∝ Frequency; Filter: Repeats ≥ 5, No NT)", fontsize=16, fontweight='bold')
    plt.axis("off")
    plt.tight_layout()
    plt.savefig("spa_MST_Strict_Sized.png", dpi=300)
    plt.show()

else:
    print(" No data suitable for analyses, empty data set")
    print("Check repeat length.")
