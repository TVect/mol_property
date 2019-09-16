# molecular similarity

ref: [分子相似性搜索的一些工具](https://www.tvect.cn/archives/394)

**Usage**

First of all, you need to build the faiss index.

1. use the following command to get molecular fingerprint vectors.

```
cd mol_property/similarity & python utils.py
```

2. use `mol_property/similarity/HammingSS.build_index` or `mol_property/similarity/CosineSS.build_index` to generate a new index file.

