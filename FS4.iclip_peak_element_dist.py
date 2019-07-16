# h4_full_12_all, count>=6
fig, ax = plt.subplots(figsize=(6,6))

sizes = [532, 4972, 16188]
labels = ['UTR5', 'CDS', 'UTR3']
total = sum(sizes)
ax.pie(sizes, labels=labels, 
        autopct=lambda(p): '{:.0f}'.format(p * total / 100), shadow=False, startangle=140)



ax.axis('equal')


plt.savefig('/Users/gongjing/SeafileSyn/Project/zebrafish/0-figures-format/iclip_tag_dist_h4_rep12.c6.pdf')
plt.show()

# h6_full_12_all, count>=6
fig, ax = plt.subplots(figsize=(6,6))

sizes = [796, 6972, 17700]
labels = ['UTR5', 'CDS', 'UTR3']
total = sum(sizes)
ax.pie(sizes, labels=labels, 
        autopct=lambda(p): '{:.0f}'.format(p * total / 100), shadow=False, startangle=140)



ax.axis('equal')


plt.savefig('/Users/gongjing/SeafileSyn/Project/zebrafish/0-figures-format/iclip_tag_dist_h6_rep12.c6.pdf')
plt.show()