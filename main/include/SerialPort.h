#pragma once
#ifndef SERIALPORT_H_
#define SERIALPORT_H_

#include <Windows.h>

/** ����ͨ����
 *   
 *  ����ʵ���˶Դ��ڵĻ�������
 *  �����������ָ�����ڵ����ݡ�����ָ�����ݵ�����
 */
class CSerialPort
{
public:
	CSerialPort(void);
	~CSerialPort(void);

public:
	
	/** ��ʼ�����ں���
	 *
	 *  @param:  UINT portNo ���ڱ��,Ĭ��ֵΪ1,��COM1,ע��,������Ҫ����9
	 *  @param:  UINT baud   ������,Ĭ��Ϊ115200,������������
	 *  @param:  char parity �Ƿ������żУ��,'Y'��ʾ��Ҫ��żУ��,'N'��ʾ����Ҫ��żУ��
	 *  @param:  UINT databits ����λ�ĸ���,Ĭ��ֵΪ8������λ
	 *  @param:  UINT stopsbits ֹͣλʹ�ø�ʽ,Ĭ��ֵΪ1
	 *  @param:  DWORD dwCommEvents Ĭ��ΪEV_RXCHAR,��ֻҪ�շ�����һ���ַ�,�����һ���¼�
	 *  @return: bool  ��ʼ���Ƿ�ɹ�
	 *  @note:   ��ʹ�����������ṩ�ĺ���ǰ,���ȵ��ñ��������д��ڵĳ�ʼ��
	 *���������� �������ṩ��һЩ���õĴ��ڲ�������,����Ҫ����������ϸ��DCB����,��ʹ�����غ���
	 *           ������������ʱ���Զ��رմ���,�������ִ�йرմ���
	 *  @see:    
	 */
	bool InitPort( UINT  portNo = 3,UINT  baud = CBR_115200, char  parity = 'E', UINT  databits = 8,
		           UINT  stopsbits = 1, DWORD dwCommEvents = EV_RXCHAR);

	/** ���ڳ�ʼ������
	 *
	 *  �������ṩֱ�Ӹ���DCB�������ô��ڲ���
	 *  @param:  UINT portNo
	 *  @param:  const LPDCB & plDCB
	 *  @return: bool  ��ʼ���Ƿ�ɹ�
	 *  @note:   �������ṩ�û��Զ���ش��ڳ�ʼ������
	 *  @see:    
	 */
	bool InitPort( UINT  portNo ,const LPDCB& plDCB );

	/** ���������߳�
	 *
	 *  �������߳���ɶԴ������ݵļ���,�������յ������ݴ�ӡ����Ļ���
	 *  @return: bool  �����Ƿ�ɹ�
	 *  @note:   ���߳��Ѿ����ڿ���״̬ʱ,����flase
	 *  @see:    
	 */
	bool OpenListenThread();

	/** �رռ����߳�
	 *
	 *  
	 *  @return: bool  �����Ƿ�ɹ�
	 *  @note:   ���ñ�������,�������ڵ��߳̽��ᱻ�ر�
	 *  @see:    
	 */
	bool CloseListenTread();

    /** �򴮿�д����
	 *
	 *  ���������е�����д�뵽����
	 *  @param:  float control ����λ
	 *  @param:  float data д������
	 *  @param:  unsigned int length ��Ҫд������ݳ���
	 *  @return: bool  �����Ƿ�ɹ�
	 *  @note:   length��Ҫ����pData��ָ�򻺳����Ĵ�С
	 *  @see:    
	 */
	bool WriteData(unsigned char control, float data, unsigned int length);

	/** ��ȡ���ڻ������е��ֽ���
	 *
	 *  
	 *  @return: UINT  �����Ƿ�ɹ�
	 *  @note:   �����ڻ�������������ʱ,����0
	 *  @see:    
	 */
	UINT GetBytesInCOM();

	/** ��ȡ���ڽ��ջ�������һ���ֽڵ�����
	 *
	 *  
	 *  @param:  char & cRecved ��Ŷ�ȡ���ݵ��ַ�����
	 *  @return: bool  ��ȡ�Ƿ�ɹ�
	 *  @note:   
	 *  @see:    
	 */
	bool ReadChar(char &cRecved);

	bool ReadByte(unsigned char* cRecved, DWORD bytesToRead);

	double GetDelta(double cmdDelta, double stime);

	bool findValidFrame(unsigned char* buffer, DWORD bufferSize, float& result);

	float readFloatData();

private:

	/** �򿪴���
	 *
	 *  
	 *  @param:  UINT portNo �����豸��
	 *  @return: bool  ���Ƿ�ɹ�
	 *  @note:   
	 *  @see:    
	 */
	bool openPort( UINT  portNo );

	/** �رմ���
	 *
	 *  
	 *  @return: void  �����Ƿ�ɹ�
	 *  @note:   
	 *  @see:    
	 */
	void ClosePort();
	
	/** ���ڼ����߳�
	 *
	 *  �������Դ��ڵ����ݺ���Ϣ
	 *  @param:  void * pParam �̲߳���
	 *  @return: UINT WINAPI �̷߳���ֵ
	 *  @note:   
	 *  @see:    
	 */
	static UINT WINAPI ListenThread(void* pParam);

private:

	/** ���ھ�� */ 
	HANDLE  m_hComm;

	/** �߳��˳���־���� */ 
	static bool s_bExit;

	/** �߳̾�� */ 
	volatile HANDLE    m_hListenThread;

	/** ͬ������,�ٽ������� */ 
	CRITICAL_SECTION   m_csCommunicationSync;       //!< �����������

	float delta;

};

#endif //SERIALPORT_H_