def harmonize_parameter_sweeping(adata, theta_param_pool, lamda_param_pool,
                                 skip_existed = True, dry_run=False, run_umap = True, **kwargs
                                ):
    """
    theta_param_pool: list of theta values to try. Larger theta values result in more diverse clusters.
    lamda_param_pool: list of lambda valeus to try. Smaller values result in more aggressive correction.
    
    parameters for harmonize() could be passed in by kwargs
    example kwargs:
    batch_key = 'batch', use_gpu=True, verbose=True,
    max_iter_harmony =100, max_iter_clustering=250, tol_clustering=1e-7, tol_harmony=1e-7,
    random_state = 19940929, 
    """
    import umap, time, numpy as np, copy
    from umap.umap_ import nearest_neighbors
    from harmony import harmonize
    assert 'X_pca' in adata.obsm.keys()
    for x in theta_param_pool: assert x>=0
    for x in lamda_param_pool: assert x>0
        
    
    adata.uns['theta_param_pool'] = theta_param_pool
    adata.uns['lamda_param_pool'] = lamda_param_pool
    for theta in theta_param_pool  :
        for lam in lamda_param_pool:
            print('\n\n-------------------------')
            embeddingZ_name  = "Z_t%s-L%s"%(theta,lam)
            embedding2d_name = "X_umap-t%s-L%s"%(theta,lam)
            if dry_run:
                print('>DRYRUN: performing integration, generating %s\t and \t%s'%(embeddingZ_name, embedding2d_name))
                continue
            
            t1 = -1;
            if (embeddingZ_name not in adata.obsm.keys()) or (not skip_existed):
                t1 = time.time(); print(">computing ", embeddingZ_name);
                Z = harmonize(adata.obsm['X_pca'], adata.obs, 
                              theta=theta, ridge_lambda=lam,
                              **kwargs)
                adata.obsm[embeddingZ_name] = copy.copy(Z)
                
            if not run_umap:
                t2 = time.time(); print(t2-t1, 'secs elapsed');
                continue
            if (embedding2d_name not in adata.obsm.keys()) or (not skip_existed):
                print(">computing ", embedding2d_name)
                knn_obj = nearest_neighbors(
                                        adata.obsm[embeddingZ_name],
                                        n_neighbors=20, metric='euclidean',
                                        metric_kwds=None, angular=False,
                                        random_state=19940929, low_memory=False,
                                        n_jobs=32
                                        )
                mapper= umap.UMAP(n_neighbors=20, min_dist    =0.15,
                                  n_jobs =32,     random_state=19940929, 
                                  densmap=False,  dens_lambda =0.1, precomputed_knn=knn_obj,
                                 ).fit( adata.obsm[embeddingZ_name] )
                x, y=mapper.embedding_[:, 0], mapper.embedding_[:, 1]
                Xumap=np.array([x, y]).T
                adata.obsm[embedding2d_name] = Xumap.copy()
                t2 = time.time()
                print(t2-t1, 'secs elapsed')
        
def clear_parameter_sweeping(adata):
    if 'theta_param_pool' in adata.uns.keys():
        del adata.uns['theta_param_pool']
    if 'lamda_param_pool' in adata.uns.keys():
        del adata.uns['lamda_param_pool']
    for k in adata.uns.keys():
        if str(k).startswith('Z_t') or str(k).startswith('X_umap-t'):
            del adata.uns[str(k)]
        
