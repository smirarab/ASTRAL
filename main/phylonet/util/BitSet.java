package phylonet.util;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.ObjectStreamField;
import java.io.Serializable;
import java.util.Arrays;

public class BitSet implements Cloneable, Serializable {

    int hash = 0;

    private static final int ADDRESS_BITS_PER_WORD = 6;
    private static final int BITS_PER_WORD = 64;
    private static final int BIT_INDEX_MASK = BITS_PER_WORD - 1;
    private static final long WORD_MASK = -1;

    // public static boolean PRINT = false;
    private static int wordIndex(int bitIndex) {
	return bitIndex >> ADDRESS_BITS_PER_WORD;
    }

    private void checkInvariants() {
	if (!$assertionsDisabled && wordsInUse != 0
		&& words[wordsInUse - 1] == 0L)
	    throw new AssertionError();
	if (!$assertionsDisabled
		&& (wordsInUse < 0 || wordsInUse > words.length))
	    throw new AssertionError();
	if (!$assertionsDisabled && wordsInUse != words.length
		&& words[wordsInUse] != 0L)
	    throw new AssertionError();
	else
	    return;
    }

    private void recalculateWordsInUse() {
	int i;
	for (i = wordsInUse - 1; i >= 0 && words[i] == 0L; i--)
	    ;
	wordsInUse = i + 1;
	hash = 0;
    }

    public BitSet() {
	wordsInUse = 0;
	sizeIsSticky = false;
	initWords(BITS_PER_WORD);
	sizeIsSticky = false;
    }

    public BitSet(int nbits) {
	wordsInUse = 0;
	sizeIsSticky = false;
	if (nbits < 0)
	    throw new NegativeArraySizeException((new StringBuilder())
		    .append("nbits < 0: ").append(nbits).toString());

	initWords(nbits);
	sizeIsSticky = true;
	return;

    }

    private void initWords(int nbits) {
	words = new long[wordIndex(nbits - 1) + 1];
    }

    private BitSet(long words[]) {
	wordsInUse = 0;
	sizeIsSticky = false;
	this.words = words;
	this.wordsInUse = words.length;
	checkInvariants();
    }

    private void ensureCapacity(int wordsRequired) {
	if (words.length < wordsRequired) {
	    int request = Math.max(2 * words.length, wordsRequired);
	    words = Arrays.copyOf(words, request);
	    sizeIsSticky = false;
	}
    }

    private void expandTo(int wordIndex) {
	int wordsRequired = wordIndex + 1;
	if (wordsInUse < wordsRequired) {
	    ensureCapacity(wordsRequired);
	    wordsInUse = wordsRequired;
	}
    }

    private static void checkRange(int fromIndex, int toIndex) {
	if (fromIndex < 0)
	    throw new IndexOutOfBoundsException((new StringBuilder())
		    .append("fromIndex < 0: ").append(fromIndex).toString());
	if (toIndex < 0)
	    throw new IndexOutOfBoundsException((new StringBuilder())
		    .append("toIndex < 0: ").append(toIndex).toString());
	if (fromIndex > toIndex)
	    throw new IndexOutOfBoundsException((new StringBuilder())
		    .append("fromIndex: ").append(fromIndex)
		    .append(" > toIndex: ").append(toIndex).toString());
	else
	    return;
    }

    public void flip(int bitIndex) {
	if (bitIndex < 0) {
	    throw new IndexOutOfBoundsException((new StringBuilder())
		    .append("bitIndex < 0: ").append(bitIndex).toString());
	} else {
	    int wordIndex = wordIndex(bitIndex);
	    expandTo(wordIndex);
	    words[wordIndex] ^= (1L << bitIndex);
	    recalculateWordsInUse();
	    checkInvariants();
	    return;
	}
    }

    public void flip(int fromIndex, int toIndex) {
	checkRange(fromIndex, toIndex);
	if (fromIndex == toIndex)
	    return;
	int startWordIndex = wordIndex(fromIndex);
	int endWordIndex = wordIndex(toIndex - 1);
	expandTo(endWordIndex);
	long firstWordMask = WORD_MASK << fromIndex;
	long lastWordMask = WORD_MASK >>> -toIndex;
	if (startWordIndex == endWordIndex) {
	    words[startWordIndex] ^= (firstWordMask & lastWordMask);
	} else {
	    words[startWordIndex] ^= firstWordMask;
	    for (int i = startWordIndex + 1; i < endWordIndex; i++)
		words[i] ^= WORD_MASK;
	    // Handle last word
	    words[endWordIndex] ^= lastWordMask;
	}
	
	recalculateWordsInUse();
	checkInvariants();
    }

    public void set(int bitIndex) {
	if (bitIndex < 0) {
	    throw new IndexOutOfBoundsException((new StringBuilder())
		    .append("bitIndex < 0: ").append(bitIndex).toString());
	} else {
	    int wordIndex = wordIndex(bitIndex);
	    expandTo(wordIndex);
	    words[wordIndex] |= (1L << bitIndex);
	    checkInvariants();
	    hash = 0;
	    return;
	}
    }

    public void set(int bitIndex, boolean value) {
	if (value)
	    set(bitIndex);
	else
	    clear(bitIndex);
	hash = 0;
    }

    public void set(int fromIndex, int toIndex) {
	checkRange(fromIndex, toIndex);
	if (fromIndex == toIndex)
	    return;
	int startWordIndex = wordIndex(fromIndex);
	int endWordIndex = wordIndex(toIndex - 1);
	expandTo(endWordIndex);
	long firstWordMask = WORD_MASK << fromIndex;
	long lastWordMask = WORD_MASK >>> -toIndex;
	if (startWordIndex == endWordIndex) {
	    words[startWordIndex] |= firstWordMask & lastWordMask;
	} else {
	    words[startWordIndex] |= firstWordMask;
	    for (int i = startWordIndex + 1; i < endWordIndex; i++)
		words[i] = WORD_MASK;
	    words[endWordIndex] |= lastWordMask;
	}
	checkInvariants();
	hash = 0;
    }

    public void set(int fromIndex, int toIndex, boolean value) {
	if (value)
	    set(fromIndex, toIndex);
	else
	    clear(fromIndex, toIndex);
	hash = 0;
    }

    public void clear(int bitIndex) {
	if (bitIndex < 0)
	    throw new IndexOutOfBoundsException((new StringBuilder())
		    .append("bitIndex < 0: ").append(bitIndex).toString());
	int wordIndex = wordIndex(bitIndex);
	if (wordIndex >= wordsInUse)
	    return;
	
	words[wordIndex] &= ~(1L << bitIndex);
	
	recalculateWordsInUse();
	checkInvariants();
    }

    public void clear(int fromIndex, int toIndex) {
	checkRange(fromIndex, toIndex);
	if (fromIndex == toIndex)
	    return;
	int startWordIndex = wordIndex(fromIndex);
	if (startWordIndex >= wordsInUse)
	    return;
	int endWordIndex = wordIndex(toIndex - 1);
	if (endWordIndex >= wordsInUse) {
	    toIndex = length();
	    endWordIndex = wordsInUse - 1;
	}
	long firstWordMask = WORD_MASK << fromIndex;
	long lastWordMask = WORD_MASK >>> -toIndex;
	if (startWordIndex == endWordIndex) {
	    words[startWordIndex] &= ~(firstWordMask & lastWordMask);
	} else {
	    words[startWordIndex] &= ~firstWordMask;
	    for (int i = startWordIndex + 1; i < endWordIndex; i++)
		words[i] = 0L;

	    words[endWordIndex] &= ~lastWordMask;
	}

	recalculateWordsInUse();
	checkInvariants();
    }

    public void clear() {
	while (wordsInUse > 0)
	    words[--wordsInUse] = 0L;
	hash = 0;
    }

    public boolean get(int bitIndex) {
	if (bitIndex < 0) {
	    throw new IndexOutOfBoundsException((new StringBuilder())
		    .append("bitIndex < 0: ").append(bitIndex).toString());
	} else {
	    checkInvariants();
	    int wordIndex = wordIndex(bitIndex);
	    return wordIndex < wordsInUse
		    && (words[wordIndex] & 1L << bitIndex) != 0L;
	}
    }

    public BitSet get(int fromIndex, int toIndex) {
	checkRange(fromIndex, toIndex);
	checkInvariants();
	int len = length();
	if (len <= fromIndex || fromIndex == toIndex)
	    return new BitSet(0);
	if (toIndex > len)
	    toIndex = len;
	BitSet result = new BitSet(toIndex - fromIndex);
	int targetWords = wordIndex(toIndex - fromIndex - 1) + 1;
	int sourceIndex = wordIndex(fromIndex);
	boolean wordAligned = (fromIndex & BIT_INDEX_MASK) == 0;
	for (int i = 0; i < targetWords - 1;) {
	    result.words[i] = wordAligned ? words[sourceIndex]: 
		(words[sourceIndex] >>> fromIndex) |
		(words[sourceIndex+1] << -fromIndex);
	    i++;
	    sourceIndex++;
	}

	long lastWordMask = WORD_MASK >>> -toIndex;
	result.words[targetWords - 1] =(toIndex-1 & BIT_INDEX_MASK) >= (fromIndex & BIT_INDEX_MASK) ? (words[sourceIndex] & lastWordMask) >>> fromIndex
		: words[sourceIndex] >>> fromIndex
			| (words[sourceIndex + 1] & lastWordMask) << -fromIndex;

	result.wordsInUse = targetWords;
	result.recalculateWordsInUse();
	result.checkInvariants();

	return result;
    }

    public int nextSetBit(int fromIndex) {
	if (fromIndex < 0)
	    throw new IndexOutOfBoundsException((new StringBuilder())
		    .append("fromIndex < 0: ").append(fromIndex).toString());
	checkInvariants();
	int u = wordIndex(fromIndex);
	if (u >= wordsInUse)
	    return -1;
	long word = words[u] & WORD_MASK << fromIndex;
	do {
	    if (word != 0L)
		return u * 64 + Long.numberOfTrailingZeros(word);
	    if (++u == wordsInUse)
		return -1;
	    word = words[u];
	} while (true);
    }

    public int nextClearBit(int fromIndex) {
	if (fromIndex < 0)
	    throw new IndexOutOfBoundsException((new StringBuilder())
		    .append("fromIndex < 0: ").append(fromIndex).toString());
	checkInvariants();
	int u = wordIndex(fromIndex);
	if (u >= wordsInUse)
	    return fromIndex;
	long word = ~words[u] & WORD_MASK << fromIndex;
	do {
	    if (word != 0L)
		return u * 64 + Long.numberOfTrailingZeros(word);
	    if (++u == wordsInUse)
		return wordsInUse * 64;
	    word = ~words[u];
	} while (true);
    }

    public int length() {
	if (wordsInUse == 0)
	    return 0;
	else
	    return 64 * (wordsInUse - 1)
		    + (64 - Long.numberOfLeadingZeros(words[wordsInUse - 1]));
    }

    public boolean isEmpty() {
	return wordsInUse == 0;
    }

    public boolean intersects(BitSet set) {
	for (int i = Math.min(wordsInUse, set.wordsInUse) - 1; i >= 0; i--)
	    if ((words[i] & set.words[i]) != 0L)
		return true;

	return false;
    }

    public int intersectionSize(BitSet set) {

	int sum = 0;
	for (int i = Math.min(wordsInUse, set.wordsInUse) - 1; i >= 0; i--)
	    sum += Long.bitCount(words[i] & set.words[i]);

	return sum;
    }

    public int cardinality() {
	int sum = 0;
	for (int i = 0; i < wordsInUse; i++)
	    sum += Long.bitCount(words[i]);
	return sum;
    }

    public void and(BitSet set) {
	if (this == set)
	    return;
	while (wordsInUse > set.wordsInUse)
	    words[--wordsInUse] = 0L;
	for (int i = 0; i < wordsInUse; i++)
	    words[i] &= set.words[i];

	recalculateWordsInUse();
	checkInvariants();
    }

    public void or(BitSet set) {
	if (this == set)
	    return;
	int wordsInCommon = Math.min(wordsInUse, set.wordsInUse);
	if (wordsInUse < set.wordsInUse) {
	    ensureCapacity(set.wordsInUse);
	    wordsInUse = set.wordsInUse;
	}
	for (int i = 0; i < wordsInCommon; i++)
	    words[i] |= set.words[i];

	if (wordsInCommon < set.wordsInUse)
	    System.arraycopy(set.words, wordsInCommon, words, wordsInCommon,
		    wordsInUse - wordsInCommon);
	hash = 0;
	checkInvariants();
    }

    public boolean contains(BitSet other) {
	if (this == other)
	    return true;
	int wordsInCommon = Math.min(wordsInUse, other.wordsInUse);
	int i;
	for (i = 0; i < wordsInCommon; i++) {
	    if ((words[i] | other.words[i]) != words[i]) {
		return false;
	    }
	}
	// if the other Bitset has extra words, they have to be zero
	for (int j = i; j < other.words.length; ++j)
	    if (other.words[j] != 0)
		return false;
	return true;
    }

    public void xor(BitSet set) {
	int wordsInCommon = Math.min(wordsInUse, set.wordsInUse);
	if (wordsInUse < set.wordsInUse) {
	    ensureCapacity(set.wordsInUse);
	    wordsInUse = set.wordsInUse;
	}
	for (int i = 0; i < wordsInCommon; i++)
	    words[i] ^= set.words[i];

	if (wordsInCommon < set.wordsInUse)
	    System.arraycopy(set.words, wordsInCommon, words, wordsInCommon,
		    set.wordsInUse - wordsInCommon);
	recalculateWordsInUse();
	checkInvariants();
    }

    public void andNot(BitSet set) {
	for (int i = Math.min(wordsInUse, set.wordsInUse) - 1; i >= 0; i--)
	    words[i] &= ~set.words[i];

	recalculateWordsInUse();
	checkInvariants();
    }

    public int hashCode() {
	if (hash == 0) {
	    long h = 1234L;
	    for (int i = wordsInUse; --i >= 0;)
		h ^= words[i] * (long) (i + 1);

	    hash = (int) (h >> 32 ^ h);
	}
	return hash;
    }

    public int size() {
	return words.length * 64;
    }

    public boolean equals(Object obj) {
	BitSet bs = (BitSet) obj;
	if (obj == this)
	    return true;
	int max = Math.min(words.length, bs.words.length);
	int i;
	for (i = 0; i < max; ++i)
	    if (words[i] != bs.words[i])
		return false;
	// If one is larger, check to make sure all extra words are 0.
	for (int j = i; j < words.length; ++j)
	    if (words[j] != 0)
		return false;
	for (int j = i; j < bs.words.length; ++j)
	    if (bs.words[j] != 0)
		return false;
	return true;
    }

    public Object clone() {
	if (!sizeIsSticky)
	    trimToSize();
	BitSet result;

	try {
	    result = (BitSet) super.clone();
	} catch (CloneNotSupportedException e) {
	    // TODO Auto-generated catch block
	    e.printStackTrace();
	    return null;
	}

	result.words = (long[]) words.clone();
	result.checkInvariants();
	result.hash = hash;
	return result;

    }

    private void trimToSize() {
	if (wordsInUse != words.length) {
	    words = Arrays.copyOf(words, wordsInUse);
	    checkInvariants();
	}
    }

    private void writeObject(ObjectOutputStream s) throws IOException {
	checkInvariants();
	if (!sizeIsSticky)
	    trimToSize();
	java.io.ObjectOutputStream.PutField fields = s.putFields();
	fields.put("bits", words);
	s.writeFields();
    }

    private void readObject(ObjectInputStream s) throws IOException,
	    ClassNotFoundException {
	java.io.ObjectInputStream.GetField fields = s.readFields();
	words = (long[]) (long[]) fields.get("bits", null);
	wordsInUse = words.length;
	recalculateWordsInUse();
	sizeIsSticky = words.length > 0 && words[words.length - 1] == 0L;
	checkInvariants();
    }

    public String toString() {
	checkInvariants();
	int numBits = wordsInUse <= 128 ? wordsInUse * 64 : cardinality();
	StringBuilder b = new StringBuilder(6 * numBits + 2);
	b.append('{');
	int i = nextSetBit(0);
	if (i != -1) {
	    b.append(i);
	    for (i = nextSetBit(i + 1); i >= 0; i = nextSetBit(i + 1)) {
		int endOfRun = nextClearBit(i);
		do
		    b.append(", ").append(i);
		while (++i < endOfRun);
	    }

	}
	b.append('}');
	return b.toString();
    }
    
    public String toString2() {
	checkInvariants();
	int numBits = wordsInUse <= 128 ? wordsInUse * 64 : cardinality();
	StringBuilder b = new StringBuilder();
	b.append('{');
	for (int i = wordsInUse - 1 ; i>=0; i--) {
		b.append(words[i]+" ");
	}
	b.append('}');
	return b.toString();
    }

    private static final ObjectStreamField serialPersistentFields[] = { new ObjectStreamField(
	    "bits", Byte.class) };
    private long words[];
    private transient int wordsInUse;
    private transient boolean sizeIsSticky;
    private static final long serialVersionUID = 7997698588986878753L;
    static final boolean $assertionsDisabled = true; // !java/util/BitSet.desiredAssertionStatus();

}

/*
 * DECOMPILATION REPORT
 * 
 * Decompiled from: /usr/lib/jvm/java-6-openjdk/jre/lib/rt.jar Total time: 106
 * ms Jad reported messages/errors: The class file version is 49.0 (only 45.3,
 * 46.0 and 47.0 are supported) Couldn't fully decompile method clone Couldn't
 * resolve all exception handlers in method clone Exit status: 0 Caught
 * exceptions:
 */
