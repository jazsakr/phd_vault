![[violin_plot_example1.png]]
```python fold title="plot code"
# Define color palette manually
color_palette = ["#5CD32C", "#2CB0D3", "#A32CD3", "#D34F2C"]

fig_title=f"{mod_name} at TSS region of genes"
path = "/path/to/the/directory/for/plots"
fig_name=f'{path}/violin.png'
fig_stats_name=f'{path}/violin.txt'
print(f"{fig_name}\n")

melted_df = merged_df.melt(
    value_vars=["Sample 1", "Sample 2", "Sample 3", "Sample 4"], 
    var_name="Sample", 
    value_name="Methylation Level"
)

plt.figure(figsize=(10, 6))

sns.violinplot(x="Sample", y="Methylation Level", data=melted_df, 
    palette=color_palette, inner="box")

# Calculate statistics for 'Methylation Level'
summary_stats = melted_df.groupby('Sample')["Methylation Level"].describe()
summary_stats = summary_stats[['25%', '50%', '75%', 'mean', 'min', 'max']]
print(f"{fig_title}\n")
print(summary_stats)

# Add titles and labels
plt.title(fig_title, fontsize=16)
plt.xlabel("Sample", fontsize=16)
plt.ylabel("Methylation Level", fontsize=16)
plt.xticks(ticks=[0,1,2,3], labels=["Sample 1", "Sample 2", "Sample 3", "Sample 4"], fontsize=14, rotation=45)
plt.yticks(fontsize=14)
plt.tight_layout()

# Save files
plt.savefig(fig_name, bbox_inches = 'tight', dpi=300)
with open(fig_stats_name, "w") as f:
    f.write(f"{fig_title}\n")
    f.write(summary_stats.to_string())

plt.show()
```

![[violin_plot_example2.png]]
```python fold title="plot code"
# Define color palette manually
color_palette = ["#5CD32C", "#2CB0D3", "#A32CD3", "#D34F2C"]

fig_title=f"{mod_name} at TSS region of genes"
path = "/path/to/the/directory/for/plots"
fig_name=f'{path}/violin.png'
fig_stats_name=f'{path}/violin.txt'
print(f"{fig_name}\n")

melted_df = merged_df.melt(
    value_vars=["Sample 1", "Sample 2", "Sample 3", "Sample 4"], 
    var_name="Sample", 
    value_name="Methylation Level"
)

plt.figure(figsize=(10, 6))

sns.violinplot(x="Sample", y="Methylation Level", data=melted_df, 
    palette=color_palette, inner="box")

# Overlay bigger box plot
sns.boxplot(x="Sample", y="Methylation Level", data=melted_df, 
    width=0.07,  # Set width of the box plot
    showcaps=False,  # Remove the caps (horizontal lines on top/bottom of the box plot)
    boxprops={'facecolor':'black', "zorder":2},  # Set box color to black and layer it properly in z-order
    medianprops={'color':'white', 'linewidth':2},  # Set the median line to white and set thickness
    showfliers=False,  # Hide outliers (individual points outside the whiskers)
    whiskerprops={'color':'black', 'linewidth':1},  # Set whiskers (vertical lines) thickness and color
    linewidth=0.5  # Set line width of the plot elements (axes, etc.)
)

# Calculate statistics for 'Methylation Level'
summary_stats = melted_df.groupby('Sample')["Methylation Level"].describe()
summary_stats = summary_stats[['25%', '50%', '75%', 'mean', 'min', 'max']]
print(f"{fig_title}\n")
print(summary_stats)

# Add titles and labels
plt.title(fig_title, fontsize=16)
plt.xlabel("Sample", fontsize=16)
plt.ylabel("Methylation Level", fontsize=16)
plt.xticks(ticks=[0,1,2,3], labels=["Sample 1", "Sample 2", "Sample 3", "Sample 4"], fontsize=14, rotation=45)
plt.yticks(fontsize=14)
plt.tight_layout()

# Save files
plt.savefig(fig_name, bbox_inches = 'tight', dpi=300)
with open(fig_stats_name, "w") as f:
    f.write(f"{fig_title}\n")
    f.write(summary_stats.to_string())

plt.show()
```