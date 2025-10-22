![[double_barplot_example1.png]]

```python fold title:"plot code"
import time
import pandas as pd
print("pandas: " + pd.__version__)
import seaborn as sns
print("seaborn: " + sns.__version__)
import matplotlib
import matplotlib.pyplot as plt
# embed plots in notebook
%matplotlib inline
print("matplotlib: " + matplotlib.__version__)
print(time.strftime("%Y-%m-%d"))

data = {
    "sample": ["control", "control", "disease", "disease", "control", "control", "disease", "disease"],
    "stage": ["condition1", "condition1", "condition1", "condition1", "condition2", "condition2", "condition2", "condition2"],
    "treatment": ["untreated", "treated", "untreated", "treated", "untreated", "treated", "untreated", "treated"],
    "count": [4937374, 4534587, 1461742, 2543375, 4781316, 4834587, 3186521, 4316045],
}
df = pd.DataFrame(data)

fig_title="Figure Title"
path = "/path/to/the/directory/for/plots"
fig_name=f'{path}/barplot.png'
print(fig_name)

# Separate data for the two subplots
df1 = df[df["stage"] == "condition1"]
df2 = df[df["stage"] == "condition2"]

# Define custom color palettes
df1_colors = {"untreated": "mediumpurple", "treated": "lightgreen"}
df2_colors = {"untreated": "indigo", "treated": "darkgreen"}

# Set up the figure and axes
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4), sharey=True)

# condition 1 plot
sns.barplot(x="sample", y="count", hue="treatment", data=df1, palette=df1_colors, ax=ax1)
ax1.set_title("Condition 1", fontsize=16)
ax1.set_xlabel("") # Remove individual x-axis label
ax1.set_ylabel("Counts ($10^6$)", fontsize=16)
ax1.tick_params(axis="x", labelsize=16)
ax1.tick_params(axis="y", labelsize=16)
#ax1.grid(False)
for bar in ax1.containers:
    ax1.bar_label(bar, fmt="%.0f", fontsize=11, label_type="edge", padding=3)

# condition 2 plot
sns.barplot(x="sample", y="count", hue="treatment", data=df2, palette=df2_colors, ax=ax2)
ax2.set_title("Condition 2", fontsize=16)
ax2.set_xlabel("") # Remove individual x-axis label
ax2.tick_params(axis="x", labelsize=16)
ax2.tick_params(axis="y", labelsize=16)
#ax2.grid(False)
for bar in ax2.containers:
    ax2.bar_label(bar, fmt="%.0f", fontsize=11, label_type="edge", padding=3)

# Get handles and labels from each plot
h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()

# Remove the default legends from the subplots
ax1.get_legend().remove()
ax2.get_legend().remove()

# Create the first legend on the right side
fig.legend(h1, l1, title="Condition 1", loc='center left', bbox_to_anchor=(0.88, 0.75), frameon=False, fontsize=14, title_fontsize=14)

# Create the second legend on the right side, below the first one
fig.legend(h2, l2, title="Condition 2", loc='center left', bbox_to_anchor=(0.88, 0.45), frameon=False, fontsize=14, title_fontsize=14)

# Add a single, centered x-axis label  and title for the entire figure
fig.supxlabel("Samples", fontsize=16, y=0.07)
fig.suptitle(fig_title, fontsize=18, y=0.95)

# Adjust layout to make space for the legends
plt.ylim(0, 5500000)
plt.tight_layout(rect=[0, 0, 0.9, 1]) # rect=[left, bottom, right, top]

# Save the plot
#plt.savefig(fig_name, bbox_inches = 'tight', dpi=300)

plt.show()
```