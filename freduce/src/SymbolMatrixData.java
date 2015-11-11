/*
 * SymbolMatrixData.java - Todd Riley
 *
 */

import java.io.*;
import java.lang.String;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.PackedSymbolList;

public class SymbolMatrixData implements Serializable {

    // original ordering of seqIDs
    String[] originalSeqIDs = null;

    ////////////////////////////////////////////////////////////////////////////////
    // seqIDs must be in same order as the data structs below!!
    ////////////////////////////////////////////////////////////////////////////////
    String[] seqIDs = null;

    // The sequences for each seqID are strored as bit-packed PackedSymbolLists
    // packedSymbolLists[seqIDNum] = PackedSymbolList
    PackedSymbolList[] packedSymbolLists = null;
}
