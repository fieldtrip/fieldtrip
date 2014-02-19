//=============================================================================================================
/**
* @file     mne_rt_buffer_reader.java
* @author   Christoph Dinh <chdinh@nmr.mgh.harvard.edu>;
*           Matti Hamalainen <msh@nmr.mgh.harvard.edu>
* @version  1.0
* @date     July, 2012
*
* @section  LICENSE
*
* Copyright (C) 2012, Christoph Dinh and Matti Hamalainen. All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that
* the following conditions are met:
*     * Redistributions of source code must retain the above copyright notice, this list of conditions and the
*       following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and
*       the following disclaimer in the documentation and/or other materials provided with the distribution.
*     * Neither the name of the Massachusetts General Hospital nor the names of its contributors may be used
*       to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
* WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
* PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL MASSACHUSETTS GENERAL HOSPITAL BE LIABLE FOR ANY DIRECT,
* INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
* HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
* NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*
*
* @brief    Contains the implementation of the MNERTBufferReader Class, to realize fast array reading.
*
*/

import java.io.*;

class mne_rt_buffer_reader
{
    public mne_rt_buffer_reader(DataInput p_DataInputStream)
    {
        m_DataInputStream = p_DataInputStream;
    }

    public byte[] readBuffer(int p_iLength)
    {
        byte[] t_byteBuffer = new byte[p_iLength];

        try
        {
            m_DataInputStream.readFully(t_byteBuffer, 0, p_iLength);
        }

        catch (StreamCorruptedException e)
        {
            System.out.println("Stream Corrupted Exception Occured");
            t_byteBuffer = new byte[0];
        }
        catch (EOFException e)
        {
            System.out.println("EOF Reached");
            t_byteBuffer = new byte[0];
        }
        catch (IOException e)
        {
            System.out.println("IO Exception Occured");
            t_byteBuffer = new byte[0];
        }

        return t_byteBuffer;
    }

    private DataInput m_DataInputStream;
}
