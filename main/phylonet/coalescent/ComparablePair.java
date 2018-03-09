package phylonet.coalescent;


public class ComparablePair<K, V extends Comparable<V>>implements Comparable{
	public K key;
	public V value;
	public ComparablePair(K key, V value) {
		this.key = key;
		this.value = value;
	}

	@Override
	public int compareTo(Object pair2) {
		return value.compareTo(((ComparablePair<K, V>)pair2).value);
	}

}
