fig, ax = plt.subplots(1,2,figsize=(6,6))

sizes = [145, 103]
labels = ['Maternal stable', 'others']
total = sum(sizes)
ax[0].pie(sizes, labels=labels, 
        autopct=lambda(p): '{:.0f}'.format(p * total / 100), shadow=False, startangle=140)
ax[0].axis('equal')


sizes = [496, 321]
labels = ['Maternal decay', 'others']
total = sum(sizes)
ax[1].pie(sizes, labels=labels, 
        autopct=lambda(p): '{:.0f}'.format(p * total / 100), shadow=False, startangle=140)
ax[1].axis('equal')

plt.savefig('/Users/gongjing/SeafileSyn/Project/zebrafish/0-figures-format/iCLIP_structure_change_materl_decay_pct.pdf')
plt.show()