# Voronoids
A parallel code to compute Voronoi diagram for future cosmological surveys

## Todo

- [x] Optimize identify non-conflict points in parallel version. Currently, we check every pair of points within the set of vertices through indexing. But I think there might be gain in recording the boundary of the sites, and just use a n-box intersection algorithm. 
- [ ] Optimize computing sphere volume and radius
- [ ] Optimize construction of neighbor relationship between the newly established points.
- [x] Pipelining threads. Currently, there is a lot of dead time between pushing data into the new tree and computing. It would be nice to separate them using channel.
- [ ] Better KDtree implementation? The KDtree package is not very performant, sometime one query on a tree with 10000 points can take 30 microseconds.
- [x] Parallel insertion? For non-conflict points, the number of vertices and simplicies are known ahead of time. So it should be possible to partition the memory and insert them in parallel. For a vector of set of conflict points, the number of simplicies being inserted probably cannot be determined ahead of time. But the respective part of the tree are independent, so maybe it can be done?
