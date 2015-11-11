/*
 * Output.java - Todd Riley
 *
 */

/*
  This class performs simple logging and message displaying to the user. It does NOT buffer the output messages!
  They are immediately sent to the destination output (at the expense of efficiency).

  1) All messages that have a level >= the minLvel are sent to the destination output.

  2) Levels have a byte value from min=0 to max=100

  3) if a message doesn't have a level then it's default value is 100.

  4) Setting minLevel=101 turns off ALL messaging!

  5) dateTimeStamp places a date-time stamp at the beginning of all the messages

    Here's a example heirarchy of message levels:

     0 = all lowest-level debug messages
     25 = low-level debug messages
     50 = medium-level debug messages
     75 = high-level debug messages
     100 = messages that should ALWAYS be logged (when messaging is turned on)

 */



import java.util.*;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.StringTokenizer;
import java.io.File;


class Output {

    enum Destination {STDOUT, STDERR, FILE, BUFFER}

    private StringBuffer buffer = null;
    private byte minLevel = 0;
    private Destination destination = null;
    private boolean dateTimeStamp = false;

    private PrintWriter printWriter = null;

    public Output(Destination aDest, byte aMinLevel, boolean aDateTimeStamp) {
        this.destination = aDest;
        this.minLevel = aMinLevel;
        this.dateTimeStamp = aDateTimeStamp;

        switch (aDest) {
        case STDOUT:
            printWriter = new PrintWriter(new OutputStreamWriter(System.out));
            break;
        case STDERR:
            printWriter = new PrintWriter(new OutputStreamWriter(System.err));
            break;
        case FILE:
            System.err.print("Error: In output.java must initialize with FilePathName when Destination is FILE.");
            break;
        case BUFFER:
            buffer = new StringBuffer();
            break;
        }
    }

    public Output(String aFilePathName, byte aMinLevel, boolean aDateTimeStamp) {
        this.destination = Destination.FILE;
        this.minLevel = aMinLevel;
        this.dateTimeStamp = aDateTimeStamp;

        try {
            printWriter = new PrintWriter(new OutputStreamWriter(new FileOutputStream(aFilePathName)));
        }
        catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void setMinLevel (byte aLevel) {
        this.minLevel = aLevel;
    }

    public void print(String aMessage) {
        this.print((byte)100, aMessage);
    }

    public void print(byte aLevel, String aMessage) {

        if (aLevel >= this.minLevel) {

            switch (this.destination) {
            case BUFFER:
                buffer.append(aMessage);
                break;
            default:
                printWriter.print(aMessage);
                printWriter.flush();
                break;
            }
        }
    }

    public String get() {
        return(buffer.toString());
    }

    public void purge() {
        buffer = new StringBuffer();
    }

    public void println(String aMessage) {
        this.print((byte)100, aMessage + "\n");
    }

    public void println(byte aLevel, String aMessage) {
        this.print(aLevel, aMessage + "\n");
    }

    public void println() {
        this.print((byte)100, "\n");
    }

    public void println(byte aLevel) {
        this.print(aLevel, "\n");
    }

}
