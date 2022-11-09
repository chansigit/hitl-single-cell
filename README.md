# Human-in-the-loop Single-cell Analysis

---
Manual selection of data is actually very important for single-cell analysis. 

Here we show how cells can be selected with interactive widgets in Jupyter.

```
import jscatter
lst = interactive_scatter_list(adata, 
                         var_list=['Cd3d','ann220413' ], 
                         embedding_list=['X_umap-t4.0-L1.0','X_umap-t4.0-L1.0'])

jscatter.link(lst, rows=1)
```

<img width="912" alt="image" src="https://user-images.githubusercontent.com/18084613/200885958-7618ee5b-feb2-4c9f-8cc3-6d320b5617f4.png">


```
# get selected cells
cell_number_1  = lst[0].selection()
cell_indices_1 = adata.obs_names[cell_number_1]
print("# of selected cells", len(cell_indices_1))
```

<img width="400" alt="image" src="https://user-images.githubusercontent.com/18084613/200886076-690fcce8-9e90-4797-a112-a39239bc6419.png">
