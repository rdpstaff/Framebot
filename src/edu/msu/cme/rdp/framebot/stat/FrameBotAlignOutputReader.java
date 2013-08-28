/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.framebot.stat;


import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

/**
 *
 * @author wangqion
 */
public class FrameBotAlignOutputReader extends FrameBotReaderCore {
    
    private Scanner scanner;
    private boolean ignoreAlignment = false ; // default only keep the STAT line, ignore the alignment
    
    public FrameBotAlignOutputReader(File framebotResult, boolean ignoreAlignment) throws FileNotFoundException{
        scanner = new Scanner(framebotResult).useDelimiter(">");
        this.ignoreAlignment = ignoreAlignment;
    }
    
    public boolean hasNext(){
        return scanner.hasNext();
    }
    
    public FrameBotStat next(){
        String[] statline = scanner.next().split("\n");
        // the second line 
        if ( ignoreAlignment) {
            return getStatLine(statline[1]);
        }else {
            return getAlignment(statline);
        }
    }
     
    public void close(){
        scanner.close();
    }
    
    
}
