
# Learning Optimal Fair Policies


This is an implementation of the following paper: 

[Learning Optimal Fair Policies](http://proceedings.mlr.press/v97/nabi19a/nabi19a.pdf)
(ICML, 2019)



If you find it useful, please consider citing:
```
@inproceedings{nabi19fairpolicy, 
title = {Learning Optimal Fair Policies},
author = { Razieh Nabi and Daniel Malinsky and Ilya Shpitser},
booktitle = {Proceedings of the Thirty Sixth International Conference on Machine Learning (ICML-36th)},
year = { 2019 }, 
publisher = {PMLR}, 
volume = {97}, 
pages = {4674-4682}
}
```


## Abstract 
Systematic discriminatory biases present in our society influence the way data is collected and stored, the way variables are defined, and the way scientific findings are put into practice as policy. Automated decision procedures and learning algorithms applied to such data may serve to perpet- uate existing injustice or unfairness in our society. In this paper, we consider how to make optimal but fair decisions, which “break the cycle of injus- tice” by correcting for the unfair dependence of both decisions and outcomes on sensitive features (e.g., variables that correspond to gender, race, disability, or other protected attributes). We use methods from causal inference and constrained optimization to learn optimal policies in a way that addresses multiple potential biases which af- flict data analysis in sensitive contexts, extending the approach of Nabi & Shpitser (2018). Our pro- posal comes equipped with the theoretical guaran- tee that the chosen fair policy will induce a joint distribution for new instances that satisfies given fairness constraints. We illustrate our approach with both synthetic data and real criminal justice data.


## License
[MIT](https://choosealicense.com/licenses/mit/)
