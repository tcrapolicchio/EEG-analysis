%prova
clear all
chanlocs = load("chanlocs.mat");
prova = load("data_example.mat");

chanlocs = chanlocs.chanlocs;
prova = prova.prova;

topoplot(prova,chanlocs,'maplimits',([0,1]));
colorbar