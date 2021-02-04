# TopicTracker

## Introduction
<p align="justify">
Topic trajectory information provides precious insight into the dynamics of topics and enlightens evolution relationships between topics over a given time. 
Nevertheless, the implementation of a large body of the existing methods for Topic trajectory identification is not available as usable software. To tackle this issue, we have developed <b>TopicTracker</b> that is a platform for topic trajectory identification and visualisation. The key of <b>TopicTracker</b> is that it can represent the three facets of information together: (a) evolutionary pathways of dynamic topics with evolution strengths, (b) evolution states of the topics and (c) the topic importance. <b>TopicTracker</b> is a publicly available software implemented as a package for the R software. 
</p>

## Setup Instruction
1. Clone the repository
```
git clone https://github.com/Yongbinkang/topicTracker.git
```
2. Install R packages: igraph, hash, plotrix, formattable
3. Install Jupyter notebook and R Kernel (IRkernel) for Jupyter (e.g., refer to [this](https://dzone.com/articles/using-r-on-jupyternbspnotebook))

## Directory structure
* The __`src/`__ directory contains the source code and demo code for <b>TopicTracker</b>:
  * The `topicTracker.R`contains the main source code for <b>TopicTracker</b> that define how to build a Topic Evolution Tree (TET) and visualise it using igraph. For more details about the code, please refer to our paper (see below).
  * The `TopicTrackerDemo.ipynb` contains the demo code on Jupyter notebook using R Kernel. Users can simply run the sequential steps using the provided data.
* The __`data/`__ directory contains two sample data used in `TopicTrackerDemo.ipynb`. 

## A TET Example
<img src="https://github.com/Yongbinkang/topicTracker/blob/main/image/tet.png" alt="Topic Evolution Tree Example" width="500" height="400">

The TET shows the follows:
 * Evolutionary pathways of 11 topics with evolution strengths, indicated by the directed edges with different edge colors, 
 * Evolution states of the topics, indicated by different node colors and sizes 
 * The topic popularity indicated by topic weight information on the y-axis.

## Citation
If you use <b>TopicTracker</b> in your research and development, please cite [TopicTracker: A Platform for Topic Trajectory Identification and Visualisation](https://arxiv.org/xxx)
```
@misc{kang2021topictracker,
      title={TopicTracker: A Platform for Topic Trajectory Identification and Visualisation}, 
      author={Yong-Bin Kang and xxx},
      year={2021},
      eprint={xxxx},
      archivePrefix={arXiv},
      primaryClass={cs.IR}
}
```
