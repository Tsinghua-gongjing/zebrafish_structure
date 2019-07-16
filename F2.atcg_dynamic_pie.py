fig, ax = plt.subplots(figsize=(6,6))

sizes = [1914704, 1561318, 1004011, 1138800]
labels = ['A', 'T', 'C', 'G']
ax.pie(sizes, labels=labels, 
        autopct='%1.2f%%', shadow=False, startangle=140)

ax.axis('equal')

plt.savefig('/Users/gongjing/SeafileSyn/Project/zebrafish/0-figures-format/all_stage_dynamic_combine_base_ratio.pie.pdf')
plt.show()