function [tree1, tree2] = averageMarginals(tree1,tree2,i,j)

avg = mean([tree1(:,i) tree2(:,j)],2);

tree1(:,i) = avg;
tree2(:,j) = avg;