/*
 * Copyright (C) 2010, Stefan Klanke
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */
 
package nl.fcdonders.fieldtrip;
 
/** Simple class for describing data blocks as used in GET_DAT and PUT_DAT requests */
public class DataDescription {
	public int nSamples;
	public int nChans;
	public int dataType;
	public int sizeBytes;
}