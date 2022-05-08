package phylonet.coalescent;

public abstract class AbstractWeightCalculatorMP<T> extends AbstractWeightCalculator<T>{


	int counter;

	public AbstractWeightCalculatorMP() {
		super(false);
	}
	
	public abstract  Long getWeight(T t);

	public int getCalculatedWeightCount() {
			return this.callcounter;
	}
	
	protected Long calculateWeight(T t) {
		return this.calculateWeight((T[])new Tripartition[] {(Tripartition) t})[0];
	}
	public abstract Long[] calculateWeight(T[] t);

	@Override
	public Object clone() throws CloneNotSupportedException {
		return super.clone();
	}
}