//#include "StdAfx.h"
#include "SerialPort.h"
#include <process.h>
#include <iostream>

/** �߳��˳���־ */ 
bool CSerialPort::s_bExit = false;
/** ������������ʱ,sleep���´β�ѯ�����ʱ��,��λ:�� */ 
const UINT SLEEP_TIME_INTERVAL = 5;

DWORD findFrameStart(unsigned char* buffer, DWORD bufferSize);

CSerialPort::CSerialPort(void)
: m_hListenThread(INVALID_HANDLE_VALUE)
{
	m_hComm = INVALID_HANDLE_VALUE;
	m_hListenThread = INVALID_HANDLE_VALUE;

	InitializeCriticalSection(&m_csCommunicationSync);
	delta = 180.0f;
}

CSerialPort::~CSerialPort(void)
{
	CloseListenTread();
	ClosePort();
	DeleteCriticalSection(&m_csCommunicationSync);
}

bool CSerialPort::InitPort( UINT portNo /*= 1*/,UINT baud /*= CBR_9600*/,char parity /*= 'N'*/,
						    UINT databits /*= 8*/, UINT stopsbits /*= 1*/,DWORD dwCommEvents /*= EV_RXCHAR*/ )
{

	/** ��ʱ����,���ƶ�����ת��Ϊ�ַ�����ʽ,�Թ���DCB�ṹ */ 
	char szDCBparam[50];
	sprintf_s(szDCBparam, "baud=%d parity=%c data=%d stop=%d", baud, parity, databits, stopsbits);

	/** ��ָ������,�ú����ڲ��Ѿ����ٽ�������,�����벻Ҫ�ӱ��� */ 
	if (!openPort(portNo))
	{
		return false;
	}

	/** �����ٽ�� */ 
	EnterCriticalSection(&m_csCommunicationSync);

	/** �Ƿ��д����� */ 
	BOOL bIsSuccess = TRUE;

    /** �ڴ˿���������������Ļ�������С,���������,��ϵͳ������Ĭ��ֵ.
	 *  �Լ����û�������Сʱ,Ҫע�������Դ�һЩ,���⻺�������
	 */
	if (bIsSuccess )
	{
		bIsSuccess = SetupComm(m_hComm,32768,32768);
	}

	/** ���ô��ڵĳ�ʱʱ��,����Ϊ0,��ʾ��ʹ�ó�ʱ���� */
	COMMTIMEOUTS  CommTimeouts;
	CommTimeouts.ReadIntervalTimeout         = 0;
	CommTimeouts.ReadTotalTimeoutMultiplier  = 0;
	CommTimeouts.ReadTotalTimeoutConstant    = 0;
	CommTimeouts.WriteTotalTimeoutMultiplier = 0;
	CommTimeouts.WriteTotalTimeoutConstant   = 0; 
	if ( bIsSuccess)
	{
		bIsSuccess = SetCommTimeouts(m_hComm, &CommTimeouts);
	}

	DCB  dcb;
	if ( bIsSuccess )
	{
		// ��ANSI�ַ���ת��ΪUNICODE�ַ���
		DWORD dwNum = MultiByteToWideChar (CP_ACP, 0, szDCBparam, -1, NULL, 0);
		wchar_t *pwText = new wchar_t[dwNum] ;
		if (!MultiByteToWideChar (CP_ACP, 0, szDCBparam, -1, pwText, dwNum))
		{
			bIsSuccess = TRUE;
		}

		/** ��ȡ��ǰ�������ò���,���ҹ��촮��DCB���� */ 
		bIsSuccess = GetCommState(m_hComm, &dcb) && BuildCommDCB(pwText, &dcb) ;
		/** ����RTS flow���� */ 
		dcb.fRtsControl = RTS_CONTROL_ENABLE; 

		/** �ͷ��ڴ�ռ� */ 
		delete [] pwText;
	}

	if ( bIsSuccess )
	{
		/** ʹ��DCB�������ô���״̬ */ 
		bIsSuccess = SetCommState(m_hComm, &dcb);
	}
		
	/**  ��մ��ڻ����� */
	PurgeComm(m_hComm, PURGE_RXCLEAR | PURGE_TXCLEAR | PURGE_RXABORT | PURGE_TXABORT);

	/** �뿪�ٽ�� */ 
	LeaveCriticalSection(&m_csCommunicationSync);

	return bIsSuccess==TRUE;
}

bool CSerialPort::InitPort( UINT portNo ,const LPDCB& plDCB )
{
	/** ��ָ������,�ú����ڲ��Ѿ����ٽ�������,�����벻Ҫ�ӱ��� */ 
	if (!openPort(portNo))
	{
		return false;
	}
	
	/** �����ٽ�� */ 
	EnterCriticalSection(&m_csCommunicationSync);

	/** ���ô��ڲ��� */ 
	if (!SetCommState(m_hComm, plDCB))
	{
		return false;
	}

	/**  ��մ��ڻ����� */
	PurgeComm(m_hComm, PURGE_RXCLEAR | PURGE_TXCLEAR | PURGE_RXABORT | PURGE_TXABORT);

	/** �뿪�ٽ�� */ 
	LeaveCriticalSection(&m_csCommunicationSync);

	return true;
}

void CSerialPort::ClosePort()
{
	/** ����д��ڱ��򿪣��ر��� */
	if( m_hComm != INVALID_HANDLE_VALUE )
	{
		CloseHandle( m_hComm );
		m_hComm = INVALID_HANDLE_VALUE;
	}
}

bool CSerialPort::openPort( UINT portNo )
{
	/** �����ٽ�� */ 
	EnterCriticalSection(&m_csCommunicationSync);

	/** �Ѵ��ڵı��ת��Ϊ�豸�� */ 
    char szPort[50];
	sprintf_s(szPort, "COM%d", portNo);

	/** ��ָ���Ĵ��� */ 
	m_hComm = CreateFileA(szPort,		                /** �豸��,COM1,COM2�� */ 
						 GENERIC_READ | GENERIC_WRITE,  /** ����ģʽ,��ͬʱ��д */   
						 0,                             /** ����ģʽ,0��ʾ������ */ 
					     NULL,							/** ��ȫ������,һ��ʹ��NULL */ 
					     OPEN_EXISTING,					/** �ò�����ʾ�豸�������,���򴴽�ʧ�� */ 
						 0,    
						 0);    

	/** �����ʧ�ܣ��ͷ���Դ������ */ 
	if (m_hComm == INVALID_HANDLE_VALUE)
	{
		std::cout << "���ڴ�ʧ��" << std::endl;
		LeaveCriticalSection(&m_csCommunicationSync);
		return false;
	}
	std::cout << "���ڴ򿪳ɹ�" << std::endl;

	/** �˳��ٽ��� */ 
	LeaveCriticalSection(&m_csCommunicationSync);

	return true;
}

bool CSerialPort::OpenListenThread()
{
	/** ����߳��Ƿ��Ѿ������� */ 
	if (m_hListenThread != INVALID_HANDLE_VALUE)
	{
		/** �߳��Ѿ����� */ 
		return false;
	}

	s_bExit = false;
	/** �߳�ID */ 
	UINT threadId;
	/** �����������ݼ����߳� */ 
	m_hListenThread = (HANDLE)_beginthreadex(NULL, 0, ListenThread, this, 0, &threadId);
	if (!m_hListenThread)
	{
		return false;
	}
	/** �����̵߳����ȼ�,������ͨ�߳� */ 
	if (!SetThreadPriority(m_hListenThread, THREAD_PRIORITY_ABOVE_NORMAL))
	{
		return false;
	}

	return true;
}

bool CSerialPort::CloseListenTread()
{	
	if (m_hListenThread != INVALID_HANDLE_VALUE)
	{
		/** ֪ͨ�߳��˳� */ 
		s_bExit = true;

		/** �ȴ��߳��˳� */ 
		Sleep(10);

		/** ���߳̾����Ч */ 
		CloseHandle( m_hListenThread );
		m_hListenThread = INVALID_HANDLE_VALUE;
	}
	return true;
}

UINT CSerialPort::GetBytesInCOM()
{
	DWORD dwError = 0;	/** ������ */ 
	COMSTAT  comstat;   /** COMSTAT�ṹ��,��¼ͨ���豸��״̬��Ϣ */ 
	memset(&comstat, 0, sizeof(COMSTAT));

	UINT BytesInQue = 0;
	/** �ڵ���ReadFile��WriteFile֮ǰ,ͨ�������������ǰ�����Ĵ����־ */ 
	if ( ClearCommError(m_hComm, &dwError, &comstat) )
	{
		BytesInQue = comstat.cbInQue; /** ��ȡ�����뻺�����е��ֽ��� */ 
	}

	return BytesInQue;
}

UINT WINAPI CSerialPort::ListenThread( void* pParam )
{
	/** �õ������ָ�� */ 
	CSerialPort *pSerialPort = reinterpret_cast<CSerialPort*>(pParam);

	// �߳�ѭ��,��ѯ��ʽ��ȡ��������
	while (!pSerialPort->s_bExit) 
	{
		UINT BytesInQue = pSerialPort->GetBytesInCOM();
		/** ����������뻺������������,����Ϣһ���ٲ�ѯ */ 
		if ( BytesInQue == 0 )
		{
			Sleep(SLEEP_TIME_INTERVAL);
			continue;
		}

		/** ��ȡ���뻺�����е����ݲ������ʾ */
		char cRecved = 0x00;
		do
		{
			cRecved = 0x00;
			if(pSerialPort->ReadChar(cRecved) == true)
			{
				//std::cout << cRecved ; 
				continue;
			}
		}while(--BytesInQue);
	}

	return 0;
}

bool CSerialPort::ReadChar( char &cRecved )
{
	BOOL  bResult     = TRUE;
	DWORD BytesRead   = 0;
	if(m_hComm == INVALID_HANDLE_VALUE)
	{
		return false;
	}

	/** �ٽ������� */ 
	EnterCriticalSection(&m_csCommunicationSync);

	/** �ӻ�������ȡһ���ֽڵ����� */ 
	bResult = ReadFile(m_hComm, &cRecved, 1, &BytesRead, NULL);
	if ((!bResult))
	{ 
		/** ��ȡ������,���Ը��ݸô�����������ԭ�� */ 
		DWORD dwError = GetLastError();

		/** ��մ��ڻ����� */ 
		PurgeComm(m_hComm, PURGE_RXCLEAR | PURGE_RXABORT);
		LeaveCriticalSection(&m_csCommunicationSync);

		return false;
	}

	/** �뿪�ٽ��� */ 
	LeaveCriticalSection(&m_csCommunicationSync);

	return (BytesRead == 1);

}

bool CSerialPort::ReadByte(unsigned char* buffer, DWORD bytesToRead) {
	BOOL bResult = TRUE;
	DWORD bytesRead = 0;

	if (m_hComm == INVALID_HANDLE_VALUE) {
		return false;
	}

	// �ٽ�������
	EnterCriticalSection(&m_csCommunicationSync);

	// �ӻ�������ȡָ�������ֽڵ�����
	bResult = ReadFile(m_hComm, buffer, bytesToRead, &bytesRead, NULL);

	/*for (DWORD i = 0; i < bytesRead; ++i) {
		printf("0x%02X ", buffer[i]);
	}
	printf("\n");*/

	if (!bResult) {
		// ��ȡ�����룬���Ը��ݸô�����������ԭ��
		DWORD dwError = GetLastError();

		std::cout << "Failed to read data, error code: " << dwError << std::endl;
		// ��մ��ڻ�����
		PurgeComm(m_hComm, PURGE_RXCLEAR | PURGE_RXABORT);
		LeaveCriticalSection(&m_csCommunicationSync);

		return false;
	}

	// �뿪�ٽ���
	LeaveCriticalSection(&m_csCommunicationSync);

	return (bytesRead == bytesToRead);
}
float convertToSingle(const unsigned char data[]) {
	// ��ȡ data �����еĵ� 4 ���� 7 ��Ԫ��
	unsigned char byteArray[4] = { data[3], data[4], data[5], data[6] };

	// �� byteArray ת��Ϊ float ����
	float* floatPtr = reinterpret_cast<float*>(byteArray);
	return *floatPtr;
}

bool CSerialPort::WriteData(unsigned char control, float data, unsigned int length )
{
	BOOL   bResult     = TRUE;
	DWORD  BytesToSend = 0;
	if(m_hComm == INVALID_HANDLE_VALUE)
	{
		std::cout << "Serial port is not opened" << std::endl;
		return false;
	}

	/** �ٽ������� */ 
	//EnterCriticalSection(&m_csCommunicationSync);

	// ����λ����֪����dsp�ĸ�����
	unsigned char stateIndex = static_cast<unsigned char>(control);

	unsigned char* pData = new unsigned char[8];
	pData[0] = 0xAA;
	pData[1] = 0x55;
	pData[2] = stateIndex;
	// ��float���͵�����ת��Ϊ�ֽ�����
	unsigned char* floatBytes = reinterpret_cast<unsigned char*>(&data);
	for (int i = 0; i < sizeof(float); ++i) {
		pData[3 + i] = floatBytes[i];
	}
	pData[7] = 0x0D;

	// ��ӡ���͵�ָ��
	float a = convertToSingle(pData);
	std::cout << (a-180) << std::endl;

	/** �򻺳���д��ָ���������� */ 
	bResult = WriteFile(m_hComm, pData, length, &BytesToSend, NULL);
	if (!bResult)  
	{
		DWORD dwError = GetLastError();
		std::cout << "Failed to write data, error code: " << dwError << std::endl;
		/** ��մ��ڻ����� */ 
		PurgeComm(m_hComm, PURGE_RXCLEAR | PURGE_RXABORT);
		//LeaveCriticalSection(&m_csCommunicationSync);

		return false;
	}

	/** �뿪�ٽ��� */ 
	//LeaveCriticalSection(&m_csCommunicationSync);

	return true;
}


double CSerialPort::GetDelta(double cmdDelta, double stime) {
	cmdDelta += 180.0;
	float fcmdDelta = (float)cmdDelta;
	WriteData(0x04, fcmdDelta, 8);
	Sleep(stime * 1000);
	unsigned char* Data = new unsigned char[8];

	// ��մ��ڽ��ջ�����
	if (!PurgeComm(m_hComm, PURGE_RXCLEAR | PURGE_RXABORT)) {
		LeaveCriticalSection(&m_csCommunicationSync);
		/*return false;*/
	}

	float res = readFloatData();
	//std::cout << res << std::endl;
	return (double)(res-180);
}

bool CSerialPort::findValidFrame(unsigned char* buffer, DWORD bufferSize, float& result) {
	// ������Ҫ8���ֽڲſ��ܰ�������֡
	if (bufferSize < 8) return false;

	// ������������֡ͷ
	for (DWORD i = 0; i <= bufferSize - 8; i++) {
		// ����Ƿ��ҵ�֡ͷ
		if (buffer[i] == 0xAA && buffer[i + 1] == 0x55 && buffer[i + 2] == 0x11) {
		//if ((buffer[i]&0xFF) == 0xAA && (buffer[i + 1]&0xFF) == 0x55 && (buffer[i + 2]&0xFF) == 0x11) {
			//if (buffer[i] == 0xAA && buffer[i + 1] == 0x55 && buffer[i + 2] == 0x11) {
			// ���֡β
			//if ((buffer[i + 7] & 0xFF) == 0x0D) {
			if (buffer[i + 7] == 0x0D) {
				for (DWORD j = i; j < 8; ++j) {
					printf("0x%02X ", buffer[j]);
				}
				printf("\n");
				// ��ȡ����������
				result = convertToSingle(&buffer[i + 3]);
				// ����true��ʾ�ҵ���Ч֡����ͨ�����÷��ؽ��
				return true;
			}
		}
	}
	return false;
}

float CSerialPort::readFloatData() {
	const int BUFFER_SIZE = 256;
	unsigned char buffer[BUFFER_SIZE] = { 0 };
	DWORD bufferPos = 0;
	float result = 0.0f;

	// ��մ��ڽ��ջ�����
	if (!PurgeComm(m_hComm, PURGE_RXCLEAR | PURGE_RXABORT)) {
		LeaveCriticalSection(&m_csCommunicationSync);
		/*return false;*/
	}

	while (true) {
		 // ���Զ�ȡ��������
		if (bufferPos < BUFFER_SIZE - 8) {
			unsigned char tempBuffer[8];
			DWORD bytesRead = 0;

			// ��ȡ���8���ֽ�
			if (ReadByte(tempBuffer, 8)) {
				// ���¶�ȡ��������ӵ�������
				memcpy(buffer + bufferPos, tempBuffer, 8);
				bufferPos += 8;
			}
		}
		//unsigned char c;
		//if (ReadByte(&c, 1)) {
		//	buffer[bufferPos++] = c;
		//}
		//// ���� buffer ���
		//if (bufferPos > BUFFER_SIZE - 8) {
		//	memmove(buffer, buffer + 1, --bufferPos);
		//}

		// ����Ƿ�����Ч֡
		if (findValidFrame(buffer, bufferPos, result)) {
			// �ҵ���Ч֡���������������Ƴ��Ѵ��������
			result = convertToSingle(buffer);
			memmove(buffer, buffer + findFrameStart(buffer, bufferPos), bufferPos - findFrameStart(buffer, bufferPos));
			bufferPos -= findFrameStart(buffer, bufferPos);
			return result;
		}

		/*if (findValidFrame(buffer, bufferPos, result)) {
			result = convertToSingle(buffer);
			DWORD start = findFrameStart(buffer, bufferPos);
			memmove(buffer, buffer + start + 8, bufferPos - (start + 8));
			bufferPos -= (start + 8);
			return result;
		}*/

		// ����������������Ƴ���ɵ�����
		if (bufferPos >= BUFFER_SIZE - 8) {
			memmove(buffer, buffer + 8, bufferPos - 8);
			bufferPos -= 8;
		}

		// ����CPUռ�ù���
		Sleep(1);
	}
}

// ���һ������е�һ����Ч֡����ʼλ��
// ���û���ҵ���Ч֡�����ػ��������ȣ���û����Ч֡��
DWORD findFrameStart(unsigned char* buffer, DWORD bufferSize) {
	for (DWORD i = 0; i <= bufferSize - 8; i++) {
		//if ((buffer[i] & 0xFF) == 0xAA && (buffer[i + 1] & 0xFF) == 0x55 && (buffer[i + 2] & 0xFF) == 0x11) {
		if (buffer[i] == 0xAA && buffer[i + 1] == 0x55 && buffer[i + 2]== 0x11) {
			//if (buffer[i] == 0xAA && buffer[i + 1] == 0x55 && buffer[i + 2] == 0x11) {
			return i;
		}
	}
	return bufferSize; // û���ҵ���Ч֡
}

