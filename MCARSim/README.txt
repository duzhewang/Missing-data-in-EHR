-Method1 folder: update beta conditioned on b_{i}

-Method2 folder: update beta not conditioned on b_{i}

-In both methods, we impute missing x using sub-model 2 directly.

———————————————————————————————————————————————————————————————————————

-folder3: update beta conditioned on b_{i}, impute the initial value of missing X_{it} from its own available observations, and update missing X_{it} using M-H algorithm. 

-folder4: update beta conditioned on b_{i}, impute the initial value of missing X_{it} by linear interpolation, and update missing X_{it} using M-H algorithm. 

-folder5: update beta NOT conditioned on b_{i}, impute the initial value of missing X_{it} by linear interpolation, and update missing X_{it} using M-H algorithm. 