#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler

# 1. Load your datasets
epitope_match_df = pd.read_csv("../Data/full_align_with_similarity.csv")  # Matching results
iedb_df = pd.read_csv("../Data/wrangled_IEDB_with_sequences.csv")    # Full IEDB dataset

# 2. Fix naming in IEDB dataframe to match epitope_match_df
iedb_df["Protein_source"] = iedb_df["Protein_source"].str.replace(
    "acetyltransferase component of pyruvate dehydrogenase complex",
    "acetyltransferase",
    regex=False
)

# 3. Compute Match_Count
match_counts = epitope_match_df["Epitope_Source"].value_counts()
epitope_match_df["Match_Count"] = epitope_match_df["Epitope_Source"].map(match_counts)

# 4. Compute Study_Count
study_counts = iedb_df["Protein_source"].value_counts()

# 5. Compute Species_Count: number of unique species matching each Epitope_Source
species_counts = epitope_match_df.groupby("Epitope_Source")["Organism_Source"].nunique()

# 6. Build epitope summary
epitope_summary = epitope_match_df.groupby("Epitope_Source").agg(
    Mean_Percent_Identity=("Percent_Identity", "mean"),
    Match_Count=("Epitope_Source", "count")
).reset_index()

# 7. Add Study_Count and Species_Count
epitope_summary["Study_Count"] = epitope_summary["Epitope_Source"].map(study_counts)
epitope_summary["Species_Count"] = epitope_summary["Epitope_Source"].map(species_counts)

# 8. Normalize match count
epitope_summary["Normalized_Match_Count"] = epitope_summary["Match_Count"] / epitope_summary["Study_Count"]

# 9. Drop any missing values
epitope_summary = epitope_summary.dropna()

# 10. Prepare feature matrix X
X = epitope_summary[["Mean_Percent_Identity", "Normalized_Match_Count", "Species_Count"]].values

#%%
# 11. Standardize features
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# 12. Fit GMM with 2 components
gmm = GaussianMixture(
    n_components=2,
    covariance_type="full",
    n_init=10,
    init_params="random",
    random_state=1
).fit(X_scaled)

# 13. Predict cluster labels
epitope_summary["Cluster"] = gmm.predict(X_scaled)

# 14. Print results
print("\nEpitope Sources and Assigned Clusters:")
print(epitope_summary[["Epitope_Source", "Cluster", "Mean_Percent_Identity", "Normalized_Match_Count", "Species_Count"]]
      .sort_values(by="Cluster")
      .to_string(index=False))

# 15. Plot 2D projection (just to visualize two axes)
plt.figure(figsize=(8,6))
plt.scatter(
    epitope_summary["Normalized_Match_Count"],
    epitope_summary["Mean_Percent_Identity"],
    c=epitope_summary["Cluster"],
    cmap="viridis",
    edgecolor="k",
    alpha=0.7
)
plt.xlabel("Normalized Match Count (Matches per Study)")
plt.ylabel("Mean Percent Identity (%)")
plt.title("GMM Clustering of Epitope Sources (Normalized + Species Diversity)")
plt.colorbar(label="Cluster Label")
plt.grid(True)
plt.show()

# 16. Bar plot: Number of species per epitope, sorted by Species_Count

# Sort epitope_summary by Species_Count
epitope_summary_sorted = epitope_summary.sort_values(by="Species_Count", ascending=False)

plt.figure(figsize=(14,6))
plt.bar(
    epitope_summary_sorted["Epitope_Source"],
    epitope_summary_sorted["Species_Count"],
    color="skyblue",
    edgecolor="black"
)
plt.xticks(rotation=90, ha="right")
plt.xlabel("Epitope Source")
plt.ylabel("Number of Unique Species Matched")
plt.title("Species Matches per Epitope Source")
plt.grid(axis='y')
plt.tight_layout()
plt.show()

# Shorten long name for clarity in the plot
epitope_match_df["Epitope_Source"] = epitope_match_df["Epitope_Source"].replace(
    "Dihydrolipoyllysine-residue acetyltransferase component of pyruvate dehydrogenase complex, mitochondrial",
    "acetyltransferase"
)

# ✅ Step 1: Compute median identity per source
median_order = (
    epitope_match_df
    .groupby("Epitope_Source")["Percent_Identity"]
    .median()
    .sort_values()
    .index
)

# ✅ Step 2: Reorder Epitope_Source by median value
epitope_match_df["Epitope_Source"] = pd.Categorical(
    epitope_match_df["Epitope_Source"],
    categories=median_order,
    ordered=True
)

# ✅ Step 3: Sort the dataframe by the new order
epitope_match_df = epitope_match_df.sort_values("Epitope_Source")

# ✅ Step 4: Plot as before
plt.figure(figsize=(10, 12))

ax = epitope_match_df.boxplot(
    column="Percent_Identity",
    by="Epitope_Source",
    vert=False,
    grid=False,
    patch_artist=True,
    boxprops=dict(facecolor="lightgreen", color="black"),
    medianprops=dict(color="red"),
    whiskerprops=dict(color="black"),
    capprops=dict(color="black"),
    flierprops=dict(markerfacecolor="red", marker="o", markersize=5, linestyle="none", markeredgecolor="black")
)

plt.yticks(fontsize=8)
plt.ylabel("Epitope Source")
plt.xlabel("Percent Identity (%)")

plt.title("")
plt.suptitle("")

# Minimalist axes
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)

ax.spines["left"].set_linewidth(1)
ax.spines["bottom"].set_linewidth(1)

plt.grid(False)
plt.tight_layout()

plt.savefig("../final plots/percent_identity_boxplot.png", dpi=300)
plt.show()
# 

# Define the condition for outliers
outlier_condition = (
    (epitope_summary["Mean_Percent_Identity"] > 70) &
    (epitope_summary["Species_Count"] == 1)
)

# Reassign cluster for outliers
epitope_summary.loc[outlier_condition, "Cluster"] = 0

# 18. Save final epitope summary to CSV
epitope_summary.to_csv("../Data/epitope_clusters.csv", index=False)

#%%
