# spa-MST-Generator

![Python Version](https://img.shields.io/badge/python-3.10.2-blue.svg)
![NetworkX](https://img.shields.io/badge/NetworkX-3.1-orange.svg)
![Matplotlib](https://img.shields.io/badge/Matplotlib-Visualization-green.svg)

A custom Python script designed to construct Minimum Spanning Trees (MST) based on the pairwise alignment of *Staphylococcus aureus* *spa* repeat sequences. This tool automates data filtering, genetic distance calculation, and advanced network visualization to elucidate clustering and population structure.

## 📌 Key Features

* **Automated Data Filtering:** Automatically excludes Non-typeable (NT) isolates and sequences with fewer than 5 repeat units to ensure clustering robustness.
* **Alignment-Based Clustering:** Calculates genetic divergence using the **Levenshtein edit distance** algorithm.
* **Advanced Network Layout:** Utilizes the Fruchterman-Reingold force-directed algorithm with custom spatial overlap-prevention to ensure clean, readable graphs.
* **Proportional Node Visualization:** Node sizes are dynamically scaled based on isolate frequency.
* **Source Tracking:** Automatically generates color-coded pie charts within nodes to represent epidemiological origins (e.g., Hospital, Dairy Farm, Slaughterhouse).

---

## 📂 Data Formatting Requirements

The script utilizes the `pandas` library to parse the input dataset (`.xlsx`). By default, the script is configured to read the data starting from the third row (`header=2`), assuming the first two rows contain preliminary metadata or title headers. 

Your dataset **must** include the following exact column names (case-sensitive) for the algorithm to properly process the parameters:

* `spa-type`: This column must contain the assigned *spa* genotype (e.g., t127, t005). Isolates that cannot be typed should be designated strictly as "NT".
* `Spa Repeats`: This column must contain the repeat succession sequence. The parsing function is designed to be highly flexible and will automatically extract individual repeat units whether they are separated by hyphens, commas, or spaces (e.g., `07-23-12-34` or `07,23,12,34`).
* `Sample Source`: This column specifies the epidemiological origin, isolation source, or host of the sample (e.g., Hospital, Dairy Farm). The script dynamically uses these categorical values to calculate proportions and generate the pie charts within the MST nodes.

---

## ⚙️ Methodology & Parameters

### Data Filtering and Preprocessing
Prior to network construction, strict exclusion criteria are applied programmatically to minimize analytical noise associated with short-sequence homoplasy. "NT" isolates and *spa* types containing fewer than five repeat units are filtered out.

### Distance Calculation
The genetic distance between any two given *spa* types is calculated using the **Levenshtein edit distance algorithm**. This algorithm determines the minimum number of operations required to transform one repeat sequence into another, accurately reflecting micro-evolutionary changes.

### MST Algorithm
A complete weighted undirected graph is first generated, where nodes represent distinct *spa* types and edges represent the calculated Levenshtein distance between them. The `minimum_spanning_tree` function from the NetworkX library is applied using **Kruskal’s algorithm** to extract the MST (the subset of edges that connects all nodes with the minimum total edge weight).

### Visualization Parameters
* **Layout:** The Fruchterman-Reingold force-directed algorithm (`spring_layout` in NetworkX) is utilized with the optimal distance parameter set to `k=1.5` and `iterations=1000` to ensure optimal node distribution. A custom overlap-prevention loop adjusts spatial coordinates based on node radii.
* **Node Sizing:** Node sizes are scaled proportionally to the square root of the isolate frequency using the formula:  
  *`radius = 0.03 + (√frequency × 0.015)`*  
  This visually represents prevalence without overwhelming the graphical space.
* **Edge Thickness:** Edge thicknesses are strictly categorized based on the exact Levenshtein distance to indicate the degree of similarity:
  * `Width 3.5`: distance ≤ 2.2
  * `Width 2.0`: 2.2 < distance ≤ 3.4
  * `Width 1.0`: 3.4 < distance ≤ 4.59
  * `Width 1.2 (Dashed)`: 4.59 < distance ≤ 5.8
  * `Width 1.2 (Dotted)`: distance > 5.8

---

## 🚀 Usage

1. Clone this repository to your local machine.
2. Ensure you have the required Python libraries installed:
   ```bash
   pip install pandas networkx matplotlib numpy openpyxl
