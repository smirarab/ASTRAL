package phylonet.coalescent;

import phylonet.tree.model.Tree;

public class Solution
{

    public Solution()
    {
    }

    public Tree getTree()
    {
        return this._st;
    }

    public Long getCoalNum()
    {
		return _totalCoals;
    }


	public void setCoalNum(Long _totalCoals) {
		this._totalCoals = _totalCoals;
	}

	public void setTree(Tree _st) {
		this._st = _st;
	}

	private Tree _st;
    private Long _totalCoals;
    int _clusterIDs[];
}