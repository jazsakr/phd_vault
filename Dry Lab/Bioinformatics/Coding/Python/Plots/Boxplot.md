![[boxplot_example1.png]]
```python fold title="plot code"
# Define color palette manually
color_palette = ["#5CD32C", "#2CB0D3", "#A32CD3", "#D34F2C"]

fig_title=f"{mod_name} at TSS region of genes"
path = "/path/to/the/directory/for/plots"
fig_name=f'{path}/boxplot.png'
fig_stats_name=f'{path}/boxplot.txt'
print(f"{fig_name}\n")

melted_df = merged_df.melt(
    value_vars=["Sample 1", "Sample 2", "Sample 3", "Sample 4"], 
    var_name="Sample", 
    value_name="Methylation Level"
)

plt.figure(figsize=(8, 6))

sns.boxplot(x="Sample", y="Methylation Level", data=melted_df, palette=color_palette,
    width=0.6,      # slightly narrower boxes
    fliersize=3     # smaller outliers
)

# calculate statistics for 'Methylation Level'
summary_stats = melted_df.groupby('Sample')["Methylation Level"].describe()
summary_stats = summary_stats[['25%', '50%', '75%', 'mean', 'min', 'max']]
print(f"{fig_title}\n")
print(summary_stats)

# Titles and labels
plt.title(fig_title, fontsize=16)
plt.xlabel("Sample", fontsize=16)
plt.ylabel("Methylation Level", fontsize=16)
plt.xticks(ticks=[0,1,2,3], labels=["Sample 1", "Sample 2", "Sample 3", "Sample 4"], fontsize=14, rotation=45)
plt.yticks(fontsize=14)
plt.tight_layout()

# save files
plt.savefig(fig_name, bbox_inches = 'tight', dpi=300)
with open(fig_stats_name, "w") as f:
    f.write(f"{fig_title}\n")
    f.write(summary_stats.to_string())

plt.show()
```
