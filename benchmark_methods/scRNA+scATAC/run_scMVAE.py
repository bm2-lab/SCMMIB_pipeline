# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.0
#   kernelspec:
#     display_name: python (benchmark)
#     language: python
#     name: benchmark
# ---

# %%
import os
import sys
import json
import glob
import time
import math

import numpy as np
import pandas as pd
import scanpy as sc
from tqdm import trange

import torch
from torch import optim
from torch.autograd import Variable
import torch.utils.data as data_utils

from sklearn.cluster import KMeans
from sklearn import metrics
from sklearn.metrics import cohen_kappa_score

from scMVAE.utilities import read_dataset, normalize, calculate_log_library_size, parameter_setting, save_checkpoint, load_checkpoint, adjust_learning_rate
from scMVAE.MVAE_model import scMVAE_Concat, scMVAE_NN, scMVAE_POE

# https://github.com/cmzuo11/scMVAE

# %%
from nvitop import Device
devices = Device.all()
memory_free = [device.memory_free() for device in devices]
gpu_index = memory_free.index(max(memory_free))
torch.cuda.set_device(gpu_index)

# %%
torch.set_num_threads(5)


# %%
def train(args, adata, adata1, model, train_index, test_index, lib_mean, lib_var, lib_mean1, lib_var1, real_groups, 
          final_rate, file_fla, Type1, Type, device, scale_factor ):

    train  = data_utils.TensorDataset( torch.from_numpy( adata.raw[train_index].X ),
                torch.from_numpy( lib_mean[train_index] ), 
                torch.from_numpy( lib_var[train_index] ),
                torch.from_numpy( lib_mean1[train_index] ), 
                torch.from_numpy( lib_var1[train_index] ),
                torch.from_numpy( adata1.raw[train_index].X ))
    train_loader  = data_utils.DataLoader( train, batch_size = args.batch_size, shuffle = True )

    test  = data_utils.TensorDataset( torch.from_numpy( adata.raw[test_index].X ),
                                              torch.from_numpy( lib_mean[test_index] ), 
                                              torch.from_numpy( lib_var[test_index] ),
                                              torch.from_numpy( lib_mean1[test_index] ), 
                                              torch.from_numpy( lib_var1[test_index] ),
                                              torch.from_numpy( adata1.raw[test_index].X ))
    test_loader   = data_utils.DataLoader( test, batch_size = len(test_index), shuffle = False )

    
    total         = data_utils.TensorDataset( torch.from_numpy( adata.raw.X  ),
                                              torch.from_numpy( adata1.raw.X ))
    total_loader  = data_utils.DataLoader( total, batch_size = args.batch_size , shuffle = False )
    
    args.max_epoch   = 500
    train_loss_list  = []

    flag_break       = 0
    epoch_count      = 0
    reco_epoch_test  = 0
    test_like_max    = 100000
    status = ""

    max_iteration = 10000
    args.epoch_per_test = 10

    params    = filter(lambda p: p.requires_grad, model.parameters())
    optimizer = optim.Adam( params, lr = args.lr, weight_decay = args.weight_decay, eps = args.eps )

    epoch     = 0
    iteration = 0
    start     = time.time()

    model.init_gmm_params( total_loader )

    with trange( args.max_epoch, disable=True ) as pbar:

        while True:

            model.train()

            epoch +=  1
            epoch_lr = adjust_learning_rate( args.lr, optimizer, epoch, final_rate, 10 )
            kl_weight = min( 1, epoch / args.anneal_epoch )

            for batch_idx, ( X1, lib_m, lib_v, lib_m1, lib_v1, X2 ) in enumerate(train_loader):

                X1, X2         = X1.float().to(device), X2.float().to(device)
                lib_m,lib_v    = lib_m.to(device),      lib_v.to(device)
                lib_m1, lib_v1 = lib_m1.to(device),     lib_v1.to(device)

                X1, X2         = Variable( X1 ),    Variable( X2 )
                lib_m, lib_v   = Variable( lib_m ), Variable( lib_v )
                lib_m1, lib_v1 = Variable( lib_m1 ),Variable( lib_v1 )

                optimizer.zero_grad()

                loss1, loss2, kl_divergence_l, kl_divergence_l1, kl_divergence_z =  \
                    model( X1.float(), X2.float(), lib_m, lib_v, lib_m1, lib_v1 )
                loss = torch.mean( ( scale_factor * loss1 + loss2 + kl_divergence_l + kl_divergence_l1) + (kl_weight*(kl_divergence_z)) )  

                loss.backward()
                optimizer.step()

                iteration += 1 

            epoch_count += 1

            if epoch % args.epoch_per_test == 0 and epoch > 0: 

                model.eval()

                with torch.no_grad():

                    for batch_idx, ( X1, lib_m, lib_v, lib_m1, lib_v1, X2 ) in enumerate(test_loader): 

                        X1, X2         = X1.float().to(device), X2.float().to(device)
                        lib_v, lib_m   = lib_v.to(device),      lib_m.to(device)
                        lib_v1, lib_m1 = lib_v1.to(device),     lib_m1.to(device)

                        X1, X2         = Variable( X1 ),     Variable( X2 )
                        lib_m, lib_v   = Variable( lib_m ),  Variable( lib_v )
                        lib_m1, lib_v1 = Variable( lib_m1 ), Variable( lib_v1 )
                

                        loss1, loss2, kl_divergence_l, kl_divergence_l1, kl_divergence_z = model( X1.float(), X2.float(), lib_m, lib_v, lib_m1, lib_v1 )
                        test_loss = torch.mean( ( scale_factor * loss1 + loss2 + kl_divergence_l + kl_divergence_l1) + (kl_weight*(kl_divergence_z)) )  

                        train_loss_list.append( test_loss.item() )

                        if math.isnan(test_loss.item()):
                            flag_break = 1
                            break

                        if test_like_max >  test_loss.item():
                            test_like_max   = test_loss.item()
                            epoch_count  = 0

                            save_checkpoint(model,fileName = '{}/model_best.pth.tar'.format(args.outdir))

                            print( str(epoch)+ "   " + str(loss.item()) +"   " + str(test_loss.item()) +"   " + 
                                   str(torch.mean(loss1).item()) +"   "+ str(torch.mean(loss2).item()) +
                                   "  kl_divergence_l:  " + str(torch.mean(kl_divergence_l).item()) + " kl_weight: " + str( kl_weight )+
                                   " kl_divergence_z: " + str( torch.mean(kl_divergence_z).item() ) )

            if epoch_count >= 30:
                reco_epoch_test = epoch
                status = " larger than 30 "
                break

            if flag_break == 1:
                reco_epoch_test = epoch
                status = " with NA "
                break

            if epoch >= args.max_epoch:
                reco_epoch_test = epoch
                status = " larger than 500 epoch "
                break
            
            if len(train_loss_list) >= 2 :
                if abs(train_loss_list[-1] - train_loss_list[-2]) / train_loss_list[-2] < 1e-4 :
                    reco_epoch_test = epoch
                    status = " training for the train dataset is converged! "
                    break

    duration = time.time() - start
    print('Finish training, total time: ' + str(duration) + 's' + " epoch: " + str(reco_epoch_test) + " status: " + status )

    load_checkpoint( '{}/model_best.pth.tar'.format(args.outdir), model, device)

    latent_z, recon_x1, norm_x1, recon_x2, norm_x2 = model.Denoise_batch(total_loader)

    if latent_z is not None:
        pd.DataFrame( latent_z, index= adata.obs_names ).to_csv( os.path.join( args.outdir, args.config["output_prefix"] + "-scMVAE-multi-latent.csv" ) ) 

    os.remove(os.path.join(output_path, "model_best.pth.tar"))

    # save imputation
    ## RNA
    imputation_rna = pd.DataFrame(recon_x1) # imputation
    imputation_rna.columns = adata.var.index.values
    imputation_rna.index = adata.obs.index.values
    imputation_rna.to_csv(os.path.join(output_path, config["output_prefix"] + "-scMVAE-imputation-rna.csv.gz"), compression='gzip')

    ## ATAC
    imputation_atac = pd.DataFrame(recon_x2) # imputation
    imputation_atac.columns = adata1.var.index.values
    imputation_atac.index = adata1.obs.index.values
    imputation_atac.to_csv(os.path.join(output_path, config["output_prefix"] + "-scMVAE-imputation-atac.csv.gz"), compression='gzip')

    # gpu_memory
    pid= os.getpid()        
    gpu_memory = pd.Series(dtype='str')

    devices = Device.all()
    for device in devices:
        processes = device.processes()    
        if pid in processes.keys():
            p=processes[pid]
            gpu_memory['device ' + str(device.index)] = p.gpu_memory_human()

    if len(gpu_memory):
        gpu_memory.to_csv(os.path.join(output_path, args.config['output_prefix'] + '-scMVAE-gpu_memory.csv'), header=["gpu_memory"])


# %%
def train_with_argas( args ):

    # Make Dir
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    
    # Load data
    adata = sc.read_h5ad(os.path.join(input_path, args.config['rna_h5ad_filename']))
    if not isinstance(adata.X, np.ndarray):
        adata.X = adata.X.toarray().astype('float32')

    adata1 = sc.read_h5ad(os.path.join(input_path, args.config['gam_h5ad_filename']))
    if not isinstance(adata1.X, np.ndarray):
        adata1.X = adata1.X.toarray().astype('float32')
    
    adata  = normalize(adata, filter_min_counts=False, size_factors=False, 
                       normalize_input=False, logtrans_input=False ) 
    adata1 = normalize(adata1, filter_min_counts=False, size_factors=False, 
                       normalize_input=False, logtrans_input=False ) 

    train_index, test_index = list(range( adata.n_obs )), list(range( adata.n_obs ))
    args.batch_size     = 64
    args.epoch_per_test = 10
    adata.obs['Group']  = ''
    lib_mean, lib_var   = calculate_log_library_size( adata.X )
    lib_mean1, lib_var1 = calculate_log_library_size( adata1.X )

    Nsample, Nfeature   = np.shape( adata.X )
    Nsample1, Nfeature1 = np.shape( adata1.X )

    device = torch.device("cuda" if args.use_cuda and torch.cuda.is_available() else "cpu")
    
    model  = scMVAE_POE ( encoder_1       = [Nfeature, 1024, 128, 128],
                      hidden_1        = 128, 
                      Z_DIMS          = 22, 
                      decoder_share   = [22, 128, 256],
                      share_hidden    = 128, 
                      decoder_1       = [128, 128, 1024], 
                      hidden_2        = 1024, 
                      encoder_l       = [ Nfeature, 128 ],
                      hidden3         = 128, 
                      encoder_2       = [Nfeature1, 1024, 128, 128], 
                      hidden_4        = 128,
                      encoder_l1      = [Nfeature1, 128], 
                      hidden3_1       = 128, 
                      decoder_2       = [128, 128, 1024],
                      hidden_5        = 1024, 
                      drop_rate       = 0.1, 
                      log_variational = True,
                      Type            = "ZINB", 
                      device          = device, 
                  n_centroids     = 22, 
                  penality        = "GMM",
                  model           = 1,  )

    args.lr           = 0.001
    args.anneal_epoch = 200

    model.to(device)
    infer_data = adata1

    train( args, adata, infer_data, model, train_index, test_index, lib_mean, lib_var, 
           lib_mean1, lib_var1, adata.obs['Group'], 0.0001, 1, "ZINB", "ZINB", device, 
           scale_factor = 4 )


# %%
input_path = sys.argv[1]
output_path = sys.argv[2]
config = json.load(open(sys.argv[3]))

# %%
parser = parameter_setting()
args   = parser.parse_args("")

args.workdir  =  input_path
args.outdir   =  output_path
args.config  =  config

# %%
train_with_argas(args)

# %%

# %%

# %%

# %%
# input_path = "/home/wsg/BM/data/test"
# output_path = "/home/wsg/BM/data/test/run_scMVAE"
# config = json.load(open("/home/wsg/BM/data/test/c1k.json"))

# %%
# parser = parameter_setting()
# args   = parser.parse_args("")

# args.workdir  =  input_path
# args.outdir   =  output_path
# args.config  =  config

# %%
# train_with_argas(args)
