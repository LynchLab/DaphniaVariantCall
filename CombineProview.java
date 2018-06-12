/*
======================================================================
CombineProview

Usage: Combine all mapgd proview files into one
 
Java -cp ./ CombineProview </Path/To/proview_files> <SampleID>
=====================================================================
Written by:                   
Xiaolong Wang @ Ocean University of China 
email: xiaolong@ouc.edu.cn
website: http://www.DNAplusPro.com
=====================================================================
In hope useful in genomics and bioinformatics studies.
This software is released under GNU/GPL license
Copyright (c) 2015, Ocean University of China
=====================================================================

*/

import java.io.*;
import java.util.Arrays;
//import java.util.*;
//import java.lang.String.*;

public class CombineProview{

	public static void main(String[] args){

try{
	
	CombineProview FAobj=new CombineProview();
	String[] ProFiles= new String[150];
	String[] Proview= new String[150]; 
	String[] rec= new String[150];
	
	String Com_Pro="";
	String DataDir=args[0];
	int L0=args[0].length();
	if (DataDir.substring(L0,L0)=="/"){DataDir=DataDir.substring(0,L0-1);}
	String SampleID=args[1];
	String output=DataDir+"/"+SampleID+"-combined.pro.txt";

	System.out.println("Reading proview file");
	int MaxNumOfProveiws=150;
	FileReader[] fr=new FileReader[150];
	BufferedReader[] br=new BufferedReader[150];
	
	//Looking for the pro.txt files in the data directory
	
	File file=new File(DataDir);
	File[] FileList = file.listFiles();
	Arrays.sort(FileList);
	String[] AllFiles=new String [FileList.length+1];
	System.out.println("Number of all files:"+FileList.length);
	  
	int Nf=0; //Number of proview files
	for (int i = 0; i < FileList.length; i++) {
	   if (FileList[i].isFile()) {
			AllFiles[i]=FileList[i].toString();
			//System.out.println(" reading File:"+i+": "+ AllFiles[i]);
			if(AllFiles[i].indexOf(".pro.txt")>0&&AllFiles[i].indexOf("combined.pro.txt")<0) {
				Nf++;//Number of proview files found
				ProFiles[Nf]=AllFiles[i];
				fr[Nf]=new FileReader(ProFiles[Nf]);
				br[Nf]=new BufferedReader(fr[Nf]);
				//System.out.println(ProFiles[Nf]+" is a proview file.");
			}else 
			{
				//System.out.println(AllFiles[i]+" is not a proview file.");
			}
	   }
	   if (FileList[i].isDirectory()) {
		System.out.println(FileList[i]+"is a sub directory, ignored.");
	   }
	}
		
	System.out.println(" Number of proview files  found: "+Nf);
	for(int j=1;j<=Nf;j++)
	{
		rec[j]=br[j].readLine();
		int P1=rec[j].indexOf("NAME"); //headline
		int P2=rec[j].indexOf("SCAFFOLDS"); //headline
		int P3=rec[j].indexOf("VERSION"); //headline
		int P4=rec[j].indexOf("FORMAT");
		int P5=rec[j].indexOf("CONCATENATED");
		
		if (P1>=0&&P2>P1&&P3>P2&&P4>P3&&P5>P4)
		{
			System.out.println("File "+j+": "+ProFiles[j]);
		}
		else
		{
			System.out.println("Error: "+ProFiles[j]+" is not a valid mapgd proview file, program will now exit. Please remove or rename this file and try again.");
			System.exit(0);
		}
	}
	
	BufferedWriter out=new BufferedWriter(new FileWriter(output));
	int i=0;
	//loop before the headline, copy the head 	
	while((rec[1]=br[1].readLine()) != null){
		
		for(int j=2;j<=Nf;j++)
		{
			rec[j]=br[j].readLine();
		}
		
		int P1=rec[1].indexOf("ID0"); //headline
		int P2=rec[1].indexOf("ID1"); //headline
		int P3=rec[1].indexOf("REF"); //headline
		int P4=rec[1].lastIndexOf("\t");
		
		if(P1>=0&&P2>P1&&P3>P2)//is headline
		{
			Com_Pro=rec[1];
			for(int j=2;j<=Nf;j++) //combine the headlines
			{
				Proview[j]=rec[j].substring(P4,rec[j].length());
				Com_Pro+=Proview[j];
			}
			System.out.println(Com_Pro+"<-- headline");
			out.write(Com_Pro+"\n");
			break;
		}
		else
		{
			System.out.println(rec[1]+"<-- head");
			out.write(rec[1]+"\n");			
		}
	}
	//loop after the headline
	while((rec[1]=br[1].readLine()) != null){
		
		for(int j=2;j<=Nf;j++)
		{
			rec[j]=br[j].readLine();
		}
	
			int Q1=rec[1].indexOf("scaffold");
			int Q2=rec[1].indexOf("/");
			int Q3=rec[1].lastIndexOf("\t");
			
			if(Q1>=0&&Q2>=10)//is data
			{
				i++;
				Com_Pro=rec[1];
				for(int j=2;j<=Nf;j++) //combine data for each line
				{
					Proview[j]=rec[j].substring(Q3,rec[j].length());
					Com_Pro+=Proview[j];
				}
				if (i%1000000==0){System.out.println(" "+i+" Lines combined.");}
				out.write(Com_Pro+"\n");
			}
			else
			{
				System.out.println(rec[1]+"--> tail");
				out.write(rec[1]+"\n");			
			}
	}
	System.out.println("    "+i+"Lines combined.");
	out.close();
		
	System.out.print("All mapgd proview files were combined and saved as:\n "+output);
		 
	}catch(Exception e){
		System.out.println("\nThis program is to combine multiple mapgd  proview output files (SampleID-*.pro.txt) into one single file named SampleID-combined.pro.txt\n"); 
		System.out.println("\nUsage: Java -cp ./ CombineProview </Path/To/proview_files> <SampleID>\n"); 
		System.out.println(e);
		e.printStackTrace();
	}
}
//public static int countOf (String s, char c) {
//    return s.length() - s.replace(c, "").length();
//}

}
