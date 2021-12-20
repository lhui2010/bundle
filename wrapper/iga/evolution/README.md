---
title:
slug:
date: 2021/12/20 tags: []
---


### Filter syntenic block with Ks
date: 2021-12-20 16:44:42
```bash
# Step1. Lift kaks to block
i=vv_at
python -m iga.evolution.ortho kaks_to_block  ${i}.ortho.kaks ${i}.anchors > ${i}.anchors.ks
# Step2. Filter block with Ks value
python -m iga.evolution.ortho select_block_by_ks --max_ks 1.2 ${i}.anchors.ks > ${i}.anchors.ks.filter
# Step3. Plot with filtered Ks
python -m jcvi.compara.synteny depth --histogram ${i}.anchors.ks.filter
python -m jcvi.graphics.dotplot ${i}.anchors.ks.filter
```