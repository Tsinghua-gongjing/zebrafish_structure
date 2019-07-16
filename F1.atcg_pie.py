fig, ax = plt.subplots(figsize=(6,6))

size = 0.4

val1 = [9251937, 8605104, 7445960, 7840986] # transcriptome
val2 = [4310065, 4020865, 3461761, 3679722] # icSHAPE

cmap = plt.get_cmap("tab20c")
outer_colors = cmap(np.arange(4)*4)
# inner_colors = cmap(np.array([1, 2, 5, 6, 9, 10]))

ax.pie(val1, radius=1, colors=outer_colors,
       wedgeprops=dict(width=size, edgecolor='w'), labels=['A', 'T', 'C', 'G'],autopct='%1.2f%%')

ax.pie(val2, radius=1-size,colors=outer_colors,
       wedgeprops=dict(width=size, edgecolor='w'),autopct='%1.2f%%')

ax.set(aspect="equal")
plt.savefig('/Users/gongjing/SeafileSyn/Project/zebrafish/0-figures-format/all_structure_base_ratio.pie_same_radius.pdf')
plt.show()