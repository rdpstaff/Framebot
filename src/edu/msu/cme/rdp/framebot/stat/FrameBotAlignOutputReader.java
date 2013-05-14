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
    
    public FrameBotAlignOutputReader(File framebotResult) throws FileNotFoundException{
        scanner = new Scanner(framebotResult).useDelimiter(">");
    }
    
    public boolean hasNext(){
        return scanner.hasNext();
    }
    
    public FrameBotStat next(){
        String[] statline = scanner.next().split("\n");
        // the second line        
        return getStatLine(statline[1]);
    }
     
    public void close(){
        scanner.close();
    }
    
    
}
