import argparse
import subprocess
import os
# 创建命令行解析器
parser = argparse.ArgumentParser(description='noimputation_mosaic_RNA_ADT')
# 添加布尔类型参数
parser.add_argument('--h5ad_path', type=str, help='h5ad_path')
parser.add_argument('--index', type=str, help='index')
parser.add_argument('--gpu_index',  type=str, help='gpu_index')

args = parser.parse_args()

h5ad_path = args.h5ad_path
index = args.index
gpu_index = int(args.gpu_index)


method = 'spglue'
target_folder = h5ad_path.replace('/home/wsg/BM/data','/home/sirm/project/SCMMIB') + '/'+method
if not os.path.isdir(target_folder):
    cmd = f'mkdir -p {target_folder}'
    subprocess.run(cmd,shell=True)
latent_embed_path =os.path.join(target_folder,f'RUN_{index}',f'{method}_latent.csv')



import os
import torch
import pandas as pd
import scanpy as sc
from SpatialGlue import SpatialGlue

from nvitop import Device
devices = Device.all()  # or Device.cuda.all()
memory_free = [device.memory_free() for device in devices]
if gpu_index =='No':
    gpu_index = memory_free.index(max(memory_free))
else:
    gpu_index = int(gpu_index)
print(gpu_index)
print(memory_free)
# Environment configuration. SpatialGlue pacakge can be implemented with either CPU or GPU. GPU acceleration is highly recommend for imporoved efficiency.
device = torch.device(f'cuda:{gpu_index}' if torch.cuda.is_available() else 'cpu')

# R.home位置
os.environ['R_HOME'] = '/usr/local/src/R-4.0.0'
from glob import glob

RNA_data_path = glob(f'{h5ad_path}/*RNA-counts.h5ad')[0]
ADT_data_path = glob(f'{h5ad_path}/*ADT-counts.h5ad')[0]
# read data
#file_fold = '/home/shaliu_fu/joint_bench/demo_notebook/input/spatial/human_lymph_node/' #please replace 'file_fold' with the download path
adata_omics1 = sc.read_h5ad(RNA_data_path)
adata_omics2 = sc.read_h5ad(ADT_data_path)
# adata_omics1 = sc.read_h5ad(f'{file_fold}/adata_RNA.h5ad')
# adata_omics2 = sc.read_h5ad(f'{file_fold}/adata_ADT.h5ad')

adata_omics1.var_names_make_unique()
adata_omics2.var_names_make_unique()
cell_name = list(adata_omics1.obs.index)
from SpatialGlue.preprocess import fix_seed
random_seed = 2022
fix_seed(random_seed)
from SpatialGlue.preprocess import preprocessing
data = preprocessing(adata_omics1, adata_omics2, datatype='SPOTS',n_neighbors=6)
model = SpatialGlue.SpatialGlue(data, datatype='SPOTS', device=device) # SPOTS和10X只有epoch的区别，一个200一个600，这里用SPOT的
# train model
output = model.train()
two_omics_latent = output['SpatialGlue'].copy()
two_omics_latent = pd.DataFrame(two_omics_latent,index=cell_name)
two_omics_latent.to_csv(latent_embed_path)
# 统计GPU
import os
pid = os.getpid()
device=devices[gpu_index]
processes = device.processes() 
p=processes[pid]
gpu_info = [p.gpu_memory_human(), p.command()]
tmp_suffix =  latent_embed_path.split('/')[-1]
gpu_info_txt = latent_embed_path.replace(tmp_suffix, 'GPU_info.txt')
with open(gpu_info_txt,'w') as f:
    for gpu_item in gpu_info:
        f.write(gpu_item+'\n')