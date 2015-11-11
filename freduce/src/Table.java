/*
 * Table.java - Todd Riley
 *
 */

import java.util.*;
import java.io.*;
import java.lang.*;
import java.util.regex.Pattern;
import java.text.RuleBasedCollator;
import org.biojava.bio.mydp.GreekCharsets;
import org.biojava.bio.mydp.ArrayTools;

public class Table implements Iterable, Serializable{

    private static final long serialVersionUID = -5884973024495428900L;

    private static String tab = "\\t";
    private static Pattern tabPattern = Pattern.compile(tab);


    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    // TO DO:
    // 1) keep state as to which column the table is currently sorted by
    // 2) if any table mod is performed, flag as unsorted
    // 3) if bindary search requested, sort table by that column if needed
    // 4) if generic add/set is called and table is sorted, then add/set while preserving sort
    // 5) add constructors with arg "sortBy"
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////
    // Globals
    ////////////////////////////////////////////////////////////////////////////////

    private List<List<Object>> tableList = new ArrayList<List<Object>>();
    private TableComparator tableComparator = new TableComparator();

    transient RuleBasedCollator ruleBasedCollator = GreekCharsets.getEnglishGreekCollator();
    int[] compareColumns = { 0 };
    int[] sortDirections = { 1 }; // ascending

    ////////////////////////////////////////////////////////////////////////////////
    // Constructors
    ////////////////////////////////////////////////////////////////////////////////

    public Table() {
    }

    // Loads the table as Strings
    public Table(String filePathName) {
        this(filePathName, "\t");
    }

    // Loads the table as Strings
    public Table(String filePathName, boolean translateGreek) {
        this(filePathName, "\t", translateGreek);
    }

    // Loads the table as Strings
    public Table(String filePathName, String delimiter) {
        loadFile(filePathName, delimiter);
    }

    // Loads the table as Strings
    public Table(String filePathName, String delimiter, boolean translateGreek) {
        loadFile(filePathName, delimiter, translateGreek);
    }

    public Table(Map aMap) {
        for (Object key : aMap.keySet()) {
            Object value = aMap.get(key);
            add(Arrays.asList(key, value));
        }
    }

    public Table(Object[] column0, Object[] column1) {
        if (column0.length != column1.length) {
            System.out.print("Error: column0 length ("+column0.length+") not equal to column1 length ("+column1.length+").");
        }
        for (int i=0; i < column0.length; i++) {
            add(Arrays.asList(column0[i], column1[i]));
        }
    }

    public Table(Object[] column0) {
        for (int i=0; i < column0.length; i++) {
            add(Arrays.asList(column0[i], new Integer(i)));
        }
    }

    public Table(Object[][] rows) {
        for (int i=0; i < rows.length; i++) {
            add(Arrays.asList(rows[i]));
        }
    }

    public Table(Vector<Vector<Object>> tableData) {
        for (int i=0; i < tableData.size(); i++) {
            ArrayList<Object> row = new ArrayList<Object>(tableData.get(i));
            add(row);
        }
    }

    // must set ruleBasedCollator since it is transient (not-serializable)
    private void readObject(ObjectInputStream inputStream)
            throws IOException, ClassNotFoundException
    {
        inputStream.defaultReadObject();
        ruleBasedCollator = GreekCharsets.getEnglishGreekCollator();
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Methods
    ////////////////////////////////////////////////////////////////////////////////

    public LinkedHashMap toMap(int keyIndex, int valueIndex) {
        LinkedHashMap map = new LinkedHashMap(tableList.size());

        for (int i=0; i < tableList.size(); i++ ) {
            List<Object> row = tableList.get(i);
            map.put(row.get(keyIndex), row.get(valueIndex));
        }
        return(map);
    }

    public class TableComparator implements Comparator<List<Object>>, Serializable {

        public int compare(List<Object> o1, List<Object> o2) {

            int compareValue = -2;

            // Loop thru the compareColumns and return the comparison if they don't match
            for (int i = 0; i < compareColumns.length; i++) {
                int compareColumn = compareColumns[i];
                int sortDirection = sortDirections[i];

                if (o1.get(compareColumn) instanceof Number) {
                    double o1Number = ((Number)o1.get(compareColumn)).doubleValue();
                    double o2Number = ((Number)o2.get(compareColumn)).doubleValue();
                    compareValue = Double.compare(o1Number, o2Number);
                    //System.out.print("\no1Number="+o1Number+" o2Number="+o2Number);
                }
                // else {
                //     compareValue = ruleBasedCollator.compare((String)o1.get(compareColumn), (String)o2.get(compareColumn));
                // }
                else if (o1.get(compareColumn) instanceof String) {
                    compareValue = ((String)o1.get(compareColumn)).compareToIgnoreCase((String)o2.get(compareColumn));
                }

                if ((compareValue != 0) || (i == (compareColumns.length-1))) {
                    return(sortDirection*compareValue);
                }
            }
            //System.out.print("\nvTableComparator Error: value="+compareValue);
            return (-2); // error
        }
    }

    public boolean isEmpty() {
        return(tableList.isEmpty());
    }

    //reverse the rows
    public void reverse() {
        List<List<Object>> newTableList = new ArrayList<List<Object>>();
        for (int i=tableList.size()-1; i >=0; i--) {
            newTableList.add(tableList.get(i));
        }
        tableList = newTableList;
    }

    // same as addRow
    public void add(List aList) {
        tableList.add(aList);
    }

    // add a row to bottom of the table
    public void addRow(List aList) {
        tableList.add(aList);
    }

    // same as addRow
    public void add(Vector aVector) {
        ArrayList<Object> aList = new ArrayList<Object>(aVector);
        tableList.add(aList);
    }

    // add a row to bottom of the table
    public void addRow(Vector aVector) {
        ArrayList<Object> aList = new ArrayList<Object>(aVector);
        tableList.add(aList);
    }

    // row index is 0-based indexed
    public void addRow(int index, List aList) {
        tableList.add(index, aList);
    }

    // add rows to bottom of the table
    public void addRows(List<List<Object>> aListOfRows) {
        tableList.addAll(aListOfRows);
    }

    // row index is 0-based indexed
    public void setRow(int index, List aList) {
        tableList.set(index, aList);
    }

    // col index is 0-based indexed
    public void addColumn(int index, List aList) {
        for (int i=0; i < aList.size(); i++) {
            tableList.get(i).add(index, aList.get(i));
        }
    }

    // col index is 0-based indexed
    public void setColumn(int index, List aList) {
        for (int i=0; i < aList.size(); i++) {
            tableList.get(i).set(index, aList.get(i));
        }
    }

    // add column to the RHS of table
    public void addColumn(List aList) {
        for (int i=0; i < aList.size(); i++) {
            tableList.get(i).add(aList.get(i));
        }
    }

    // ascending order column-sort
    // compareColumn is 0-based
    public void sort(int... aCompareColumns) { // aCompareColumns is really an array!!
        // int[] aCompareColumns = new int[1];
        // aCompareColumns[0] = aCompareColumn;

        sort(aCompareColumns, true); // aCompareColumns is really an array!!
    }

    // compareColumn is 0-based
    public void sort(int aCompareColumn, boolean ascending) {
        int[] aCompareColumns = new int[1];
        aCompareColumns[0] = aCompareColumn;

        sort(aCompareColumns, ascending);
    }

    // compareColumn is 0-based
    public void sort(int[] aCompareColumns, boolean ascending) {
        boolean[] ascendingArray = new boolean[aCompareColumns.length];
        for (int i=0; i<ascendingArray.length; i++) {
            ascendingArray[i] = ascending;
        }

        sort(aCompareColumns, ascendingArray);
    }

    public void sort(int[] aCompareColumns, boolean[] ascendingArray) {
        int[] sortDirectionArray = new int[aCompareColumns.length];
        for (int i=0; i<sortDirectionArray.length; i++) {
            if (ascendingArray[i]) {
                sortDirectionArray[i] = 1;
            }
            else {
                sortDirectionArray[i] = -1;
            }
        }

        sort(aCompareColumns, sortDirectionArray);
    }

    public void sort(int[] aCompareColumns, int[] directionsArray) {
        this.compareColumns = aCompareColumns;
        this.sortDirections = directionsArray;

        Collections.sort(tableList, tableComparator);
    }

    // column-based binary search requires that an ascending sort ON THAT COLUMN has been performed first!!!
    // compareColumn is 0-based
    // negative value means that the key was not found
    // if not found then value = (-(insertion point) - 1)
    //
    // NOTE : If more than one row contains the key, then ANY of those rows (with the same key) may be returned!!!
    //        If you want the first row that contains this key, then call getFirstOccurrence()
    //
    public int binarySearch(int aCompareColumn, Object key) {
        int[] compareColumns = {aCompareColumn};
        Object[] keys = {key};
        return(binarySearch(compareColumns, keys));
    }

    // column-based binary search requires that an ascending sort ON THAT COLUMN has been performed first!!!
    // compareColumn is 0-based
    // negative value means that the key was not found
    // if not found then value = (-(insertion point) - 1)
    //
    // NOTE : If more than one row contains the key, then ANY of those rows (with the same key) may be returned!!!
    //        If you want the first row that contains this key, then call getFirstOccurrence()
    //
    public int binarySearch(int[] aCompareColumns, Object[] keys) {

        if (rows() == 0) {
            return(-1);
        }

        // compareColumns = new int[1];
        // compareColumn[0] = aCompareColumn;

        // Get a copy of the first row
        List<Object> firstRow = (List<Object>)tableList.get(0) ;
        ArrayList<Object> firstRowClone = new ArrayList<Object>(firstRow);
//         ArrayList<Object> firstRow = (ArrayList<Object>)tableList.get(0) ;
//         ArrayList<Object> firstRowClone = (ArrayList<Object>)firstRow.clone();

        // set the compareColumn to the key
        for (int i=0; i < aCompareColumns.length; i++) {
            firstRowClone.set(aCompareColumns[i], keys[i]);
        }

        return(Collections.binarySearch(tableList, firstRowClone, tableComparator));

    }

    // column-based less than requires that an ascending sort ON THAT COLUMN has been performed first!!!
    // compareColumn is 0-based
    // can't use binarySearch since "If the list contains multiple elements equal to the specified object, there is no guarantee which one will be found."
    // Currently aCompareColumn must be a Number object
    public int lessThanCount(int aCompareColumn, double value) {
        for (int i=0; i < tableList.size(); i++) {
            if (((Number)tableList.get(i).get(aCompareColumn)).doubleValue() >= value) {
                return(i);
            }
        }
        return(tableList.size());
    }

    // row number is 0-based indexing
    public void removeRow(int row) {
        tableList.remove(row);
    }

    // 0-based indexing
    // begin inclusive, end exclusive
    public void removeRows(int begin, int end) {
        //((ArrayList)tableList).removeRange(begin, end);
        for (int row = end-1; row >= begin; row--) {
            removeRow(row);
        }
    }

    // 0-based indexing
    // The removal Table must be in the same ordering as Table and contain the same objects with certain rows removed!!!
    public void removeRows(Table aRemovalTable) {
        int removalTableIndex = aRemovalTable.rows() -1;

        for (int thisTableIndex = rows()-1; thisTableIndex >= 0; thisTableIndex--) {
            // List thisRow = getRow(thisTableIndex);
            List removalRow = aRemovalTable.getRow(removalTableIndex);

            if (rowEquals(thisTableIndex, removalRow)) {
                removeRow(thisTableIndex);
                removalTableIndex--;
            }
            else {
                // do nothing!
            }

            if (removalTableIndex < 0) {
                // we're done!
                break;
            }
        }
    }

    // 0-based indexing
    public boolean rowEquals(int rowIndex, List aRow) {
        List<Object> ourRow = getRow(rowIndex);
        return(ourRow.equals(aRow)); // element-wise equals
    }

    // 0-based indexing
    public boolean columnEquals(int columnIndex, List aColumn) {
        List<Object> ourColumn = getColumn(columnIndex);
        return(ourColumn.equals(aColumn)); // element-wise equals
    }

    // column number is 0-based indexing
    public void removeColumn(int column) {
        for (int i=0; i < tableList.size(); i++) {
            tableList.get(i).remove(column);
        }
    }

    // 0-based indexing
    // begin inclusive, end exclusive
    public void removeColumns(int begin, int end) {
        //((ArrayList)tableList).removeRange(begin, end);
        for (int column = end-1; column >= begin; column--) {
            removeColumn(column);
        }
    }

    // row number is 0-based indexed
    public List<Object> getRow(int rowNum) {
        return(tableList.get(rowNum));
    }

    // row, col number is 0-based indexed
    public Object getElement(int rowNum, int colNum) {
        return(tableList.get(rowNum).get(colNum));
    }

    // row, col number is 0-based indexed
    public void setElement(int rowNum, int colNum, Object element) {
        tableList.get(rowNum).set(colNum, element);
    }

    // 0-based indexing
    public List<Object> getColumn(int colNum) {
        return(getColumn(colNum, 0, tableList.size()));
    }

    // 0-based indexing
    // Performs a sort
    public List<Object> getSortedColumn(int colNum, boolean ascending) {
        sort(colNum, ascending);
        return(getColumn(colNum, 0, tableList.size()));
    }

    // 0-based indexing
    // Returns only the unique entries in the sorted column
    // Performs a sort
    public List<Object> getSortedColumnUnique(int colNum, boolean ascending) {
        sort(colNum, ascending);
        return(getColumnUnique(colNum, 0, tableList.size()));
    }

    // col number is 0-based indexed
    // begin inclusive, end exclusive
    public List<Object> getColumn(int colNum, int begin, int end) {
        List<Object> columnList = new ArrayList<Object>();
        for (int i=begin; i < end; i++) {
            columnList.add(tableList.get(i).get(colNum));
        }
        return(columnList);
    }

    // col number is 0-based indexed
    // getFirstOccurrence() requires that a sort ON THAT COLUMN has been performed first!!!
    public int getFirstOccurrence(int colNum, Object key) {
        int[] colNums = {colNum};
        Object[] keys = {key};
        return(getFirstOccurrence(colNums, keys));
    }

    public int getFirstOccurrence(int[] colNums, Object[] keys) {
        int binarySearchIndex = binarySearch(colNums, keys);

        // if not found
        if (binarySearchIndex < 0) {
            return(binarySearchIndex);
        }

        // Now traverse backwards until (colValue != key)
        for (int i=binarySearchIndex; i > 0; i--) {

            // traverse thru all the comparison columns looking for a mismtach
            for (int j=0; j < colNums.length; j++) {
                Object object = tableList.get(i).get(colNums[j]);

                if (keys[j] instanceof Number) {
                    if (((Number)object).doubleValue() != ((Number)keys[j]).doubleValue()) {
                        return(i+1);
                    }
                }
                else if (keys[j] instanceof String) {
                    if ( !((String)object).equals((String)keys[j]) ) {
                        return(i+1);
                    }
                }
            }
        }
        return(0);
    }


    // col number is 0-based indexed
    // getMatches() requires that a sort ON THAT COLUMN has been performed first!!!
    // returns null if there are no matches
    public List<List<Object>> getMatches(int colNum, Object key) {
        int[] colNums = {colNum};
        Object[] keys = {key};
        return(getMatches(colNums, keys));
    }

    // col number is 0-based indexed
    // getMatches() requires that a sort ON THAT COLUMN has been performed first!!!
    // returns null if there are no matches
    public List<List<Object>> getMatches(int[] colNums, Object[] keys) {

        List<List<Object>> matchesList = new ArrayList<List<Object>>();

        // find row of first match to search criteria
        int firstMatchRow = getFirstOccurrence(
            colNums,
            keys);


        if (firstMatchRow < 0) { // not found!
            return(null);
        }
        else { // found at least one match!

            matchesList = new ArrayList<List<Object>>();

            // get how many rows match search criteria
            int matchCount = getMatchCount(
                colNums,
                keys,
                firstMatchRow);

            // System.out.println("\nfirstMatchRow="+firstMatchRow);
            // System.out.println("matchCount="+matchCount);

            // Grab rows from the Master Table
            for (int rowIncrement = 0; rowIncrement < matchCount; rowIncrement++) {
                List rowList = getRow(firstMatchRow + rowIncrement);
                // matchesList.add(Arrays.copyOf(rowList, rowList.size()));
                matchesList.add(new ArrayList(rowList));
            }

        }

        return(matchesList);

    }

    // col number is 0-based indexed
    // getMatchCount() requires that a sort ON THAT COLUMN has been performed first!!!
    public int getMatchCount(int colNum, Object key) {
        int[] colNums = {colNum};
        Object[] keys = {key};
        return(getMatchCount(colNums, keys, 0, tableList.size()));
    }

    public int getMatchCount(int[] colNums, Object[] keys) {
        return(getMatchCount(colNums, keys, 0, tableList.size()));
    }

    // col number is 0-based indexed
    // getMatchCount() requires that a sort ON THAT COLUMN has been performed first!!!
    public int getMatchCount(int colNum, Object key, int begin) {
        int[] colNums = {colNum};
        Object[] keys = {key};
        return(getMatchCount(colNums, keys, begin, tableList.size()));
    }

    public int getMatchCount(int[] colNums, Object[] keys, int begin) {
        return(getMatchCount(colNums, keys, begin, tableList.size()));
    }

    // col number is 0-based indexed
    // getMatchCount() requires that a sort ON THAT COLUMN has been performed first!!!
    public int getMatchCount(int colNum, Object key, int begin, int end) {
        int[] colNums = {colNum};
        Object[] keys = {key};
        return(getMatchCount(colNums, keys, begin, tableList.size()));
    }

    // col number is 0-based indexed
    // getMatchCount() requires that a sort ON THAT COLUMN has been performed first!!!
    public int getMatchCount(int[] colNums, Object[] keys, int begin, int end) {
        int count = 0;
        boolean[] foundArray = new boolean[keys.length]; // have we found this key yet?
        boolean found = false;

        for (int i=begin; i < end; i++) {

            // traverse thru all the comparison columns looking for a mismtach
            for (int j=0; j < colNums.length; j++) {
                Object object = tableList.get(i).get(colNums[j]);

                if (keys[j] instanceof Number) {
                    if (((Number)object).doubleValue() == ((Number)keys[j]).doubleValue()) {
                        foundArray[j] = true;
                    }
                    else if (found) { // we've traversed passed the identical entries
                        return(count);
                    }
                    else {
                        foundArray[j] = false;
                    }
                }
                else if (keys[j] instanceof String) {
                    if ( ((String)object).equals((String)keys[j]) ) {
                        foundArray[j] = true;
                    }
                    else if (found) { // we've traversed passed the identical entries
                        return(count);
                    }
                    else {
                        foundArray[j] = false;
                    }
                }
            }

            // if we've gotten this far then maybe all keys have been simultaneuously found
            if (!ArrayTools.contains(foundArray, false)) {
                found = true;
                count++;
            }

        }
        return(count);
    }

    // getColumnUnique() requires that a sort ON THAT COLUMN has been performed first!!!
    // col number is 0-based indexed
    // Returns only the unique entries in the column
    public List<Object> getColumnUnique(int colNum) {
        return(getColumnUnique(colNum, 0, tableList.size()));
    }

    // getColumnUnique() requires that a sort ON THAT COLUMN has been performed first!!!
    // col number is 0-based indexed
    // begin inclusive, end exclusive
    // Returns only the unique entries in the column
    public List<Object> getColumnUnique(int colNum, int begin, int end) {
        List<Object> uniqueList = new ArrayList<Object>();
        Object lastObject = null;

        for (int i=begin; i < end; i++) {

            Object object = tableList.get(i).get(colNum);

            if (lastObject == null) {
                uniqueList.add(object);
                lastObject = object;
            }
            else if (object instanceof Number) {
                if (((Number)object).doubleValue() != ((Number)lastObject).doubleValue()) {
                    uniqueList.add(object);
                    lastObject = object;
                }
            }
            else if (object instanceof String) {
                if ( !((String)object).equals((String)lastObject) ) {
                    uniqueList.add(object);
                    lastObject = object;
                }
            }
        }
        return(uniqueList);
    }

    public Vector<Vector<Object>> getDataVector() {
        Vector<Vector<Object>> tableData = new Vector<Vector<Object>>();
        for (int i=0; i < rows(); i++) {
            Vector rowVector = new Vector(getRow(i));
            tableData.add(rowVector);
        }
        return(tableData);
    }

    // iterates over each row
    public Iterator<List<Object>> iterator() {
        return(tableList.iterator());
    }

    public int rows() {
        return(tableList.size());
    }

    // assumes that the table is not jagged (all rows are the same length)
    public int columns() {
        if (rows() == 0) {
            return(0);
        }
        return(tableList.get(0).size());
    }

    // get all rows
    public String toString() {
        return(toString(0, this.rows(), "\t", false));
    }

    // get all rows
    public String toString(String delimiter) {
        return(toString(0, this.rows(), delimiter, false));
    }

    // get rows in (startRow, endRow) inclusive and 0-based
    public String toString(int startRow, int endRow, String delimiter) {
        return(toString(startRow, endRow, delimiter, false));
    }

    // get rows in (startRow, endRow) inclusive and 0-based
    public String toString(int startRow, int endRow, String delimiter, boolean addRowCount) {
        return(toString(startRow, endRow, delimiter, addRowCount, false, -1));
    }

    public void scaleColumnBy(int column, double scalar) {
        scaleColumnBy(column, 0, this.rows(), scalar);
    }

    public void scaleColumnBy(int column, int startRow, int endRow, double scalar) {

        if (startRow < 0 || startRow > endRow) {
            System.err.println("Invalid Indexes: startRow="+startRow+" endRow="+endRow+".");
            return;
        }

        if (rows() == 0) {
            return;
        }

        for (int i=startRow; i < endRow; i++ ) {

            if (i < rows()) {

                List<Object> row = tableList.get(i);

                    double value = -1;

                    Object valueObject = row.get(column);
                    if (valueObject instanceof Long) {
                        value = ((Long)valueObject).doubleValue();
                    }
                    else if (valueObject instanceof Double) {
                        value = ((Double)valueObject).doubleValue();
                    }
                    else if (valueObject instanceof Float) {
                        value = ((Float)valueObject).doubleValue();
                    }
                    else if (valueObject instanceof Short) {
                        value = ((Short)valueObject).doubleValue();
                    }

                    row.set(column, new Double(value * scalar));
            }
        }
    }

    public String toString(
        int startRow,
        int endRow,
        String delimiter,
        boolean addRowCount,
        boolean createEmptyEntries,
        int createEmptyUntilExcl)
    {

        //System.out.println("Hi");

        if (startRow < 0 || startRow > endRow) {
            System.err.println("Invalid Indexes: startRow="+startRow+" endRow="+endRow+".");
            return(null);
        }

        if (rows() == 0) {
            return(null);
        }

        StringBuffer output = new StringBuffer();
        String[] lastNonEmptyEntries = new String[columns()];

        //for (int i=0; i < 1; i++ ) {
        for (int i=startRow; i < endRow; i++ ) {

            if (i < rows()) {

                List<Object> row = tableList.get(i);

                if (addRowCount) {
                    output.append((i+1)+""+delimiter);
                }

                for (int j=0; j < row.size(); j++) {

                    if (j > 0) {
                        output.append(delimiter);
                    }

                    String appendString = null;
                    if (createEmptyEntries && (j < createEmptyUntilExcl)) {
                        // If a field is identical to the last non-empty entry then replace with null
                        if (lastNonEmptyEntries[j] == null) {
                            lastNonEmptyEntries[j] = new String((String)row.get(j).toString());
                            appendString = (String)row.get(j).toString();
                        }
                        else {
                            if (((String)row.get(j).toString()).equals((String)lastNonEmptyEntries[j])) {
                                appendString = "";
                            }
                            else {
                                //lastNonEmptyEntries[j] = (String)row.get(j);
                                lastNonEmptyEntries[j] = new String((String)row.get(j).toString());
                                appendString = (String)row.get(j).toString();
                            }
                        }
                    }
                    else {
                        if (row.get(j) == null) {
                            appendString = "null";
                        }
                        else {
                            appendString = (String)row.get(j).toString();
                        }
                    }

                    output.append(appendString);
                }
                output.append("\n");
            }
        }
        return(output.toString());
    }

    // Loads the table as Strings
    public void loadFile(String filePathName) {
        loadFile(filePathName, this.tab, false);
    }

    // Loads the table as Strings
    public void loadFile(String filePathName, String delimiter) {
        loadFile(filePathName, delimiter, false);
    }

    // Loads the table as Strings
    public void loadFile(String filePathName, String delimiter, boolean translateGreek) {

        try {
            String aLine;
            BufferedReader inBuffer = new BufferedReader(new InputStreamReader(new FileInputStream(filePathName)));
            String commentToken = null;
            Pattern pattern = Pattern.compile(delimiter);

            System.out.print("\nPopulating a Table with the entries found in "+filePathName+"...");

            String[] lastNonEmptyEntries = null;

            while(((aLine = inBuffer.readLine()) != null)) {
                // include this line if line doesn't start with commentToken
                if ((commentToken == null) || !aLine.trim().startsWith(commentToken)) {

                    if (translateGreek) {
                        aLine = GreekCharsets.replaceLongGreekText(aLine);
                    }

                    String lineEntries[] = pattern.split(aLine);

                    // Fill in empty line entries with last non-empty entries
                    if (lastNonEmptyEntries == null) {
                        lastNonEmptyEntries = Arrays.copyOf(lineEntries, lineEntries.length);
                    }
                    else {
                        for (int i=0; i < lineEntries.length; i++) {
                            if ((lineEntries[i] == null) || (lineEntries[i].equals(""))) {
                                lineEntries[i] = lastNonEmptyEntries[i];
                            }
                            else {
                                lastNonEmptyEntries[i] = lineEntries[i];
                            }
                        }
                    }


                    add(Arrays.asList(lineEntries));
                }
            }
            inBuffer.close();
            System.out.println( "Done.");
        }
        catch (IOException e) {
            e.printStackTrace();
        }

    }



//     public String toString(int start, int end, String delimiter) {
//         StringBuffer output = new StringBuffer();
//         for (List<Object> row : tableList) {
//             for (Object entry : row) {
//                 if (output.length() > 0) {
//                     output.append(delimiter);
//                 }
//                 output.append(entry);
//             }
//             output.append("\n");
//         }
//         return(new String(output));
//     }

}
