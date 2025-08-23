# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 09:43:56 2024

@author: 张子慧
"""

import tifffile
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




def get_pca_feature(sample_dir,r):
    mkdir(sample_dir+'input') 
    mkdir(sample_dir+'resnet50') 
    spot= pd.read_csv(sample_dir+'/exp_location.txt',sep='\\s+',header=None)
    spot_name= pd.read_csv(sample_dir+'/spot.txt',sep='\\s+',header=None)

    #### prepare data ####
    print("---  prepare data  ---")
    centers = np.loadtxt(sample_dir+"exp_location.txt")
    img = cv2.imread(sample_dir+"tissue_hires_image.png")
    for i in range(centers.shape[0]):
        #print(i)
        x, y = int(centers[i,0]), int(centers[i,1])##### x y!
        center_name=spot_name.iloc[i,0]
        tile = get_local_image(img, x=x, y=y, r=r)
        #print(tile)
        save_dir = sample_dir+"input/{name}_r={size}.npy".format(name=center_name,size=r)
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
        
    img_list = os.listdir(sample_dir+"input/")
    #img_list = os.path.basename(x) for x in glob.glob(sample_dir+"input/*r=48.npy")]
    data1 = img_dataset1(img_dir=sample_dir+"input/", img_list=img_list,
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

    img_list = os.listdir(sample_dir+"input/")
    data1 = img_dataset(img_dir=sample_dir+"input/", img_list=img_list,
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
    np.savetxt(sample_dir+"resnet50/features.txt", features_array)
    np.savetxt(sample_dir+"resnet50/feature_name.txt", feature_name_array, fmt='%s')

    ###calculate the feature
    pca = PCA(n_components=10)
    pca.fit(features_array)
    features_pca = pca.transform(features_array)
    np.savetxt(sample_dir+"resnet50/features_pca.txt", features_pca)
    


