package phylonet.coalescent;

import phylonet.tree.model.Tree;

public class Solution
{

    public Solution()
    {
    }

    public Tree getTree()
    {
        return _st;
    }

    public Long getCoalNum()
    {
        return _totalCoals;
    }

    Tree _st;
    Long _totalCoals;
    int _clusterIDs[];
}