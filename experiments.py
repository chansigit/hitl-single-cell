import pandas as pd,numpy as np,scanpy as sc
import jscatter,scipy

def interactive_scatter_list(adata, var_list=[], embedding_list=['X_umap'], 
                       ):
    gene_list = pd.Series([v for v in var_list if v in adata.var_names])
    mat = adata[:,gene_list].X
    import scipy
    if scipy.sparse.issparse(mat): mat= mat.todense();
    expr_df   = pd.DataFrame(mat, 
                             index = adata.obs_names, 
                             columns="_expr_"+gene_list)
    
    list_of_df = []
    for emb_name in embedding_list:
        df = pd.DataFrame(adata.obsm[emb_name][:,:2],columns=[emb_name+".1", emb_name+".2"])
        list_of_df.append(df)
    emb_df = pd.concat(list_of_df, axis=1)
    emb_df.set_index(adata.obs_names, inplace=True)
    

    df = pd.concat([expr_df, emb_df, adata.obs],axis=1)
    df.set_index(adata.obs_names)
    
    list_of_scatters = []
    for var_name, emb_name in zip(var_list, embedding_list):
        sca=jscatter.Scatter(data=df, x=emb_name+'.1', y=emb_name+'.2')
        if var_name in adata.var_names:
            var_name = "_expr_"+var_name
            sca.color(by=var_name) #TODO: color configs vmin,vmax,palette 
        else:
            sca.color(by=var_name) #TODO: color configs palette 
        list_of_scatters.append(sca)
    
    return list_of_scatters
    
lst = interactive_scatter_list(adata, 
                         var_list=['Cd3d'], 
                         embedding_list=['X_umap-t4.0-L1.0'])
                         
jscatter.link(lst, rows=1)
