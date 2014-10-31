package phylonet.coalescent;

public class WQClusterCollection extends AbstractClusterCollection{

	public WQClusterCollection(int len) {
		initialize(len);
	}

	@Override
	public AbstractClusterCollection newInstance(int size) {
		return new WQClusterCollection(size);
	}
		

}
