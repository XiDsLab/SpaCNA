# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 09:43:56 2024

@author: zzh
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import time
import glob

import torch
import torchvision
import torchvision.transforms as transforms
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader

from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import cv2
import pandas as pd 

from sklearn.neighbors import NearestNeighbors
from scipy.stats import pearsonr

def mkdir(path):
 
	folder = os.path.exists(path)
 
	if not folder:                   
		os.makedirs(path)           
		print ("---  new folder...  ---")
		print ("---  OK  ---")
 
	else:
		print("---  There is this folder!  ---") 

def get_local_image(img, x, y, r):
    xmin, ymin = 0, 0
    xmax, ymax = img.shape[0], img.shape[1]
    x, y = int(x), int(y)
    if x+r>xmax or y+r>ymax or x-r<xmin or y-r<ymin:
        return np.zeros(shape=(2*r,2*r,img.shape[2]))
    return img[int(x-r):int(x+r),int(y-r):int(y+r),:]


def draw_graph(adj, centers, img, N, title, edge_color='black', save_path=None):
    plt.figure(figsize=(10, 10))
    plt.imshow(img)
    # plt.scatter(centers[:, 1], centers[:, 0], c='skyblue', s=10)
    for i in range(N):
        for j in range(i + 1, N):
            if adj[i, j] == 1:
                y0, x0 = centers[i]
                y1, x1 = centers[j]
                plt.plot([x0, x1], [y0, y1], c=edge_color, linewidth=1.5)
    plt.scatter(centers[:, 1], centers[:, 0], c='skyblue', s=10, zorder=2)

    plt.title(title)
    plt.axis('off')
    # saving
    if save_path is not None:
        plt.savefig(save_path, bbox_inches='tight', dpi=300, format='pdf')
    plt.show()



def get_pca_feature(sample_dir,r,n_neighbors = 6):
    mkdir(sample_dir+'/input') 
    mkdir(sample_dir+'/resnet50') 
    mkdir(sample_dir+'/plot') 
    #### prepare data ####
    print("---  prepare data  ---")
    # spot= pd.read_csv(sample_dir+'/exp_location.txt',sep='\\s+',header=None)
    spot_name= pd.read_csv(sample_dir+'/spot.txt',sep='\\s+',header=None)
    centers = np.loadtxt(sample_dir+"/exp_location.txt")
    img = cv2.imread(sample_dir+"/tissue_hires_image.png")
    for i in range(centers.shape[0]):
        #print(i)
        x, y = int(centers[i,0]), int(centers[i,1])##### x y!
        center_name=spot_name.iloc[i,0]
        tile = get_local_image(img, x=x, y=y, r=r)
        #print(tile)
        save_dir = sample_dir+"/input/{name}_r={size}.npy".format(name=center_name,size=r)
        np.save(save_dir, tile)
    
    #### calculate normalization mean & std ####
    img_transform1 = torchvision.transforms.Compose([
        transforms.ToPILImage(),
        transforms.Resize(224),
        transforms.CenterCrop(224),
        transforms.ToTensor()
    ])

    class img_dataset1(Dataset):
        def __init__(self, img_dir, img_list, transform=img_transform1):
            self.img_dir = img_dir
            self.img_list = img_list
            self.transform = transform

        def __getitem__(self, index):
            img_name = self.img_list[index]
            img_path = self.img_dir + img_name
            x = np.load(img_path)       
            x = self.transform(x)
            return x

        def __len__(self):
            return len(self.img_list)
        
    img_list = os.listdir(sample_dir+"/input/")
    #img_list = os.path.basename(x) for x in glob.glob(sample_dir+"input/*r=48.npy")]
    data1 = img_dataset1(img_dir=sample_dir+"/input/", img_list=img_list,
                            transform=img_transform1)
    # test = OL_data[100]
    dataloader1 = DataLoader(data1, batch_size=128, shuffle=False, num_workers=0)

    mean = 0.
    std = 0.
    nb_samples = 0.
    i = 0 
    
    for data in dataloader1:
        i += 1
        if i > 100:
            break
        batch_samples = data.size(0)
        data = data.view(batch_samples, data.size(1), -1)
        mean += data.mean(2).sum(0)
        std += data.std(2).sum(0)
        nb_samples += batch_samples

    mean /= nb_samples
    std /= nb_samples
    print(f"mean: {mean}")
    print(f"std: {std}")
    
    #### construct dataloader ####
    sample_mean = mean
    sample_std = std

    img_transform = torchvision.transforms.Compose([
        transforms.ToPILImage(),
        # transforms.Resize(256),
        transforms.CenterCrop(224),
        transforms.ToTensor(),
        transforms.Normalize(mean=sample_mean, std=sample_std),
    ])

    class img_dataset(Dataset):
        def __init__(self, img_dir, img_list, transform=img_transform):
            self.img_dir = img_dir
            self.img_list = img_list
            self.transform = transform

        def get_location(self, img_name):
            tmp = img_name.split(".")[0]
            tmp = tmp.split("_")
            name = tmp[0]
    #         pixel_x = int(tmp[0])
    #         pixel_y = int(tmp[1])
            pixel_range = int(tmp[1].split("=")[1])
            return name, pixel_range

        def __getitem__(self, index):
            img_name = self.img_list[index]
            name, _ = self.get_location(img_name)

            img_path = self.img_dir + img_name
            x = np.load(img_path)       
            x = self.transform(x)
            return x, name


        def __len__(self):
            return len(self.img_list)

    img_list = os.listdir(sample_dir+"/input/")
    data1 = img_dataset(img_dir=sample_dir+"/input/", img_list=img_list,
                            transform=img_transform)
    # test = OL_data[100]
    dataloader = DataLoader(data1, batch_size=128, shuffle=False, num_workers=0)
    
    
    #### run feature extraction model ####
    print("--- run feature extraction model resnet50  ---")
    resnet50 = torchvision.models.resnet50(pretrained=True)
    model1 = nn.Sequential(*list(resnet50.children())[:-1])
    model1.eval()
    
    features = []
    feature_name = []
    with torch.no_grad():
        for j, batch_data in enumerate(dataloader):
            #input_x, location_x, location_y = batch_data
            input_x, name = batch_data
            output_f = model1(input_x)
            output_f = torch.squeeze(output_f).detach().numpy()
            features.append(output_f)
            feature_name.append(np.array(name))

    ####save the result
    print("---save the result  ---")
    features_array = np.vstack(features)
    feature_name_array = np.hstack(feature_name).T
    np.savetxt(sample_dir+"/resnet50/features.txt", features_array)
    np.savetxt(sample_dir+"/resnet50/feature_name.txt", feature_name_array, fmt='%s')

    ###calculate the feature
    pca = PCA(n_components=10)
    pca.fit(features_array)
    features_pca = pca.transform(features_array)
    np.savetxt(sample_dir+"/resnet50/features_pca.txt", features_pca)

    #### plot ####
    df = pd.DataFrame(features_pca)
    df.index=feature_name_array
    df_sorted = df.loc[np.array(spot_name).flatten()]
    sorted_features_pca = df_sorted.values 
    sorted_feature_name_array = np.array(spot_name) 

    nbrs = NearestNeighbors(n_neighbors=n_neighbors + 1).fit(centers)
    distances, indices = nbrs.kneighbors(centers)

    N = centers.shape[0]
    distance_threshold=3*r
    adj_spatial = np.zeros((N, N), dtype=int)
    for i in range(N):
        for j in range(1, n_neighbors + 1):
            neighbor = indices[i, j]
            if distances[i, j] <= distance_threshold:
                adj_spatial[i, neighbor] = 1
                adj_spatial[neighbor, i] = 1

    corr_thresholds = [-1, -0.1, 0.1, 0.3, 0.5]
    for corr_threshold in corr_thresholds:
        adj_filtered = np.zeros((N, N), dtype=int)
        for i in range(N):
            for j in range(i + 1, N):
                if adj_spatial[i, j] == 1:
                    corr, _ = pearsonr(sorted_features_pca[i], sorted_features_pca[j])
                    if corr > corr_threshold:
                        adj_filtered[i, j] = 1
                        adj_filtered[j, i] = 1
    
        title = f"Spatial Neighbors + Feature Correlation > {corr_threshold}"
        save_path = f'{sample_dir}/plot/filtered_graph_{corr_threshold}.pdf'
        
        draw_graph(
            adj_filtered, 
            centers, 
            img, 
            N,
            title, 
            edge_color='black',
            save_path=save_path
        )
        
        print(f"processing: corr_threshold = {corr_threshold} -> {save_path}")
        

import argparse
import os
import time
from pathlib import Path

def main():
    """
    Example run:
    Suppose your data is stored in /your sample dir/,
    which contains:
        - tissue_hires_image.png   High-resolution tissue image.
        - spot.txt                 A text file with one column, listing cell barcodes (one barcode per line).
        - exp_location.txt         A text file with two columns (x, y), 
                                   representing the coordinates of each spot on the tissue_hires_image.
                                   Example:
                                       1184.9774666   1452.28040696
                                       1349.70135544   575.04009492
                                       522.94917032   1368.35209208
                                       1055.44227228  469.61972024
    
    The script will automatically generate:
        - sample1/input/            stores local tile .npy files
        - sample1/resnet50/         stores extracted features (features.txt / features_pca.txt, etc.)
    """
    parser = argparse.ArgumentParser(description='Extract spatial features from tissue images')
    parser.add_argument('--sample_dir', type=str, required=True,
                       help='Sample directory containing tissue_hires_image.png, spot.txt, exp_location.txt')
    parser.add_argument('--radius', '-r', type=int, default=50,
                       help='Cropping radius for each spot (default: 50)')
    
    args = parser.parse_args()
    
    # Use pathlib for cross-platform path handling
    sample_dir = Path(args.sample_dir)
    r = args.radius
    
    start_time = time.time()
    get_pca_feature(str(sample_dir), r)
    end_time = time.time()
    
    print(f"✅ Feature extraction completed, elapsed time: {end_time - start_time:.2f} seconds")
    # Use pathlib for path joining
    results_path = sample_dir / "resnet50"
    print(f"Results saved in: {results_path}")

if __name__ == "__main__":
    main()

