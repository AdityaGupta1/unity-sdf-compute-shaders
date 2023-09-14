#define FILL_SDFS_ON_GPU
#define FILL_SDFS_MULTIPLE_TIMES

using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using UnityEngine;
using Random = UnityEngine.Random;

enum Blocks
{
    AIR, STONE, STICK, CANDY
}

enum SdfType
{
    SPHERE, LOLLIPOP, MUSHROOM
}

struct SdfInfo
{
    public int type;
    public Vector3Int pos;
}

[RequireComponent(typeof(MeshFilter), typeof(MeshRenderer))]
public class VoxelBuilder : MonoBehaviour
{
    public ComputeShader voxelCompute;
    private int kernelHandle;

    public bool useCustomShader = false;
    public Shader meshShader;

    public bool reset = false;

    private ComputeBuffer voxelBuffer;
    private ComputeBuffer sdfInfoBuffer;

    private const int xSize = 256;
    private const int ySize = 256;
    private const int zSize = 256;

    private const int numVoxels = xSize * ySize * zSize;
    private const int numSdfs = 64;
    private const int maxNumSdfs = 256;

#if FILL_SDFS_MULTIPLE_TIMES
    private const int multipleFillSkip = 5;
    private const int multipleFillActual = 20;
#endif

    int[] blocks = new int[numVoxels];
    SdfInfo[] sdfInfos = new SdfInfo[numSdfs];

    private Mesh voxelMesh;

    private void Start()
    {
        if (voxelBuffer != null)
        {
            voxelBuffer.Release();
        }

        // cost of these functions (until GenerateSdfInfos()) should be amortized since they only happen once
        kernelHandle = voxelCompute.FindKernel("GenerateVoxels");

        voxelBuffer = new ComputeBuffer(numVoxels, sizeof(int));
        voxelCompute.SetBuffer(kernelHandle, "blocks", voxelBuffer);

        // this maxNumSdfs setup (having one permanent device array for storing SdfInfos rather than allocating a new array for each zone)
        // may not be necessary but it definitely reduces the number of memory allocations, which is probably a good thing even if it comes
        // at the cost of a bit of extra memory being allocated
        sdfInfoBuffer = new ComputeBuffer(maxNumSdfs, Marshal.SizeOf(typeof(SdfInfo)));
        voxelCompute.SetBuffer(kernelHandle, "sdfInfos", sdfInfoBuffer);

        GenerateSdfInfos();

#if FILL_SDFS_MULTIPLE_TIMES
        int totalMs = 0;
        for (int i = 0; i < multipleFillSkip + multipleFillActual; ++i)
        {
            Array.Clear(blocks, 0, numVoxels);
            int ms = FillSdfs();
            if (i >= multipleFillSkip)
            {
                totalMs += ms;
            }
        }
        Debug.LogFormat("<color=aqua>average time to fill SDFs: {0} ms</color>", totalMs / (float)multipleFillActual);
#else
        FillSdfs();
#endif

        GenerateMesh();
    }

    private void Update()
    {
        if (reset)
        {
            reset = false;
            Array.Clear(blocks, 0, numVoxels);
            int ms = FillSdfs();
            GenerateMesh();
            Debug.LogFormat("<color=green>time to fill SDFs: {0} ms</color>", ms);
        }
    }

    private void GenerateSdfInfos()
    {
        Random.InitState(76);

        for (int i = 0; i < numSdfs; ++i)
        {
            //SdfType sdfType = (i % 2 == 0) ? SdfType.SPHERE : SdfType.LOLLIPOP;
            SdfType sdfType = SdfType.MUSHROOM;
            sdfInfos[i].type = (int)sdfType;
            sdfInfos[i].pos = new Vector3Int(Random.Range(0, xSize), 0, Random.Range(0, zSize));
        }
    }

    private int FillSdfs()
    {
        DateTime t1 = DateTime.Now;

#if FILL_SDFS_ON_GPU
        sdfInfoBuffer.SetData(sdfInfos, 0, 0, numSdfs);
        voxelCompute.SetInt("numSdfs", numSdfs);

        voxelCompute.SetInt("xSize", xSize);
        voxelCompute.SetInt("ySize", ySize);
        voxelCompute.SetInt("zSize", zSize);

        voxelCompute.Dispatch(kernelHandle, xSize / 16, ySize / 16, zSize / 4);

        //voxelBuffer.GetData(blocks, 0, 0, numVoxels);
        voxelBuffer.GetData(blocks, 0, 0, xSize * 85 * zSize);
#else
        FillSdfsRecursively();
#endif

        DateTime t2 = DateTime.Now;
        int ms = t2.Subtract(t1).Milliseconds;
        Debug.LogFormat("time taken to fill SDFs: {0} ms", ms);
        return ms;
    }

#if !FILL_SDFS_ON_GPU
    float sdSphere(Vector3 p, float s) {
        return p.magnitude - s;
    }

    float sdCappedCylinder(Vector3 p, float h, float r)
    {
        Vector2 j = new Vector2(Mathf.Sqrt(p.x * p.x + p.z * p.z), p.y);
        Vector2 d = new Vector2(Mathf.Abs(j.x) - r, Mathf.Abs(j.y) - h);
        return Mathf.Min(Mathf.Max(d.x, d.y), 0) + Vector2.Max(d, Vector2.zero).magnitude;
    }

    float opUnion(float d1, float d2)
    {
        return Mathf.Min(d1, d2);
    }

    private void FillSdfsRecursively()
    {
        for (int iSdf = 0; iSdf < numSdfs; ++iSdf)
        {
            var sdfInfo = sdfInfos[iSdf];

            var start = sdfInfo.pos;
            SdfType sdfType = (SdfType)sdfInfo.type;

            HashSet<Vector3Int> visited = new HashSet<Vector3Int>();
            Queue<Vector3Int> frontier = new Queue<Vector3Int>();

            frontier.Enqueue(Vector3Int.zero);
            visited.Add(Vector3Int.zero);

            float height = 33 + Random.value * 40;

            Vector3 startPoint = Vector3.zero;
            Vector3 endPoint = startPoint + new Vector3(0, height, 0);
            Vector3[] ctrlPts = new Vector3[5];
            for (int i = 0; i < 5; ++i)
            {
                ctrlPts[i] = Vector3.Lerp(startPoint, endPoint, ((float)i) / 4);
                if (i != 0)
                {
                    Vector3 temp = new Vector3(Random.value, Random.value, Random.value) * 2 - new Vector3(1, 1, 1);
                    temp.x *= 6;
                    temp.y *= 2;
                    temp.z *= 6;
                    ctrlPts[i] += temp;
                }
            }

            Vector3[] spline = new Vector3[5];
            for (int i = 0; i < 5; ++i)
            {
                var ctrlPtsCopy = new Vector3[5];
                Array.Copy(ctrlPts, ctrlPtsCopy, 5);
                int points = 5;

                float t = ((float)i) / 4;

                while (points > 1)
                {
                    for (int j = 0; j < points - 1; ++j)
                    {
                        ctrlPtsCopy[j] = Vector3.Lerp(ctrlPtsCopy[j], ctrlPtsCopy[j + 1], t);
                    }

                    --points;
                }

                spline[i] = ctrlPtsCopy[0];
            }

            while (frontier.Count != 0)
            {
                Blocks block = Blocks.AIR;

                var posLocal = frontier.Dequeue();
                var posWorld = start + posLocal;

                visited.Add(posLocal);

                if (posWorld.x < 0 || posWorld.x >= xSize || posWorld.y < 0 || posWorld.y >= ySize || posWorld.z < 0 || posWorld.z >= zSize)
                {
                    continue;
                }

                Vector3 pos = posLocal;

                if (pos.y < -1 || pos.y > height + 12 || (Mathf.Sqrt(pos.x * pos.x + pos.z * pos.z) > 8 && (pos.y < height - 12 || (pos - new Vector3(0, height, 0)).magnitude > 35)))
                {
                    continue;
                }

                for (int i = 0; i < 5; ++i)
                {
                    var pos1 = spline[i];
                    Vector3 pos2;

                    if (i < 4)
                    {
                        pos2 = spline[i + 1];

                        if (pos.y < pos1.y - 3 || pos.y > pos2.y + 3)
                        {
                            continue;
                        }
                    }
                    else
                    {
                        pos2 = pos1 + (pos1 - spline[i - 1]).normalized * (2.5f + /*Random.value*/0.5f * 0.5f);
                    }

                    var vecLine = pos2 - pos1;

                    var pointPos = pos - pos1;
                    float ratio = Vector3.Dot(pointPos, vecLine) / Vector3.Dot(vecLine, vecLine);

                    var pointLine = vecLine * ratio;
                    var vecPointPos = pointPos - pointLine;

                    float radius;
                    Blocks potentialBlockType;
                    if (i < 4)
                    {
                        float t = (i + Mathf.Min(1, Mathf.Max(0, ratio))) / 4;
                        float x = t - 0.5f;
                        radius = 4f * x * x + 1.5f;
                        potentialBlockType = Blocks.STICK;
                    }
                    else
                    {
                        radius = (7f * /*random.y*/0.5f + 14f) * (0.8f + 0.4f * ((height - 33f) / 40f));
                        potentialBlockType = Blocks.CANDY;
                    }

                    if ((ratio >= 0 && ratio <= 1 && vecPointPos.magnitude <= radius) || (i < 4 && ratio < 0 && Vector3.Distance(pos, pos1) < radius) || (i < 3 && ratio > 1 && Vector3.Distance(pos, pos2) < radius))
                    {
                        block = potentialBlockType;
                        break;
                    }
                }

                if (block != Blocks.AIR)
                {
                    blocks[PosToIndex(posWorld.x, posWorld.y, posWorld.z)] = (int)block;

                    foreach (var direction in directions)
                    {
                        var newPos = posLocal + direction;

                        if (!visited.Contains(newPos))
                        {
                            frontier.Enqueue(newPos);
                            visited.Add(newPos);
                        }
                    }
                }
            }
        }
    }
#endif

    private int PosToIndex(int x, int y, int z)
    {
        return x + z * xSize + y * xSize * zSize;
    }

    private readonly Vector3Int[] directions =
    {
        Vector3Int.forward,
        Vector3Int.back,
        Vector3Int.right,
        Vector3Int.left,
        Vector3Int.up,
        Vector3Int.down
    };

    private readonly Vector3[][] directionVerts =
    {
        new Vector3[] { new Vector3(-0.5f, -0.5f, 0.5f), new Vector3(0.5f, -0.5f, 0.5f), new Vector3(0.5f, 0.5f, 0.5f), new Vector3(-0.5f, 0.5f, 0.5f) },
        new Vector3[] { new Vector3(0.5f, -0.5f, -0.5f), new Vector3(-0.5f, -0.5f, -0.5f), new Vector3(-0.5f, 0.5f, -0.5f), new Vector3(0.5f, 0.5f, -0.5f) },
        new Vector3[] { new Vector3(0.5f, -0.5f, 0.5f), new Vector3(0.5f, -0.5f, -0.5f), new Vector3(0.5f, 0.5f, -0.5f), new Vector3(0.5f, 0.5f, 0.5f) },
        new Vector3[] { new Vector3(-0.5f, -0.5f, -0.5f), new Vector3(-0.5f, -0.5f, 0.5f), new Vector3(-0.5f, 0.5f, 0.5f), new Vector3(-0.5f, 0.5f, -0.5f) },
        new Vector3[] { new Vector3(-0.5f, 0.5f, 0.5f), new Vector3(0.5f, 0.5f, 0.5f), new Vector3(0.5f, 0.5f, -0.5f), new Vector3(-0.5f, 0.5f, -0.5f) },
        new Vector3[] { new Vector3(-0.5f, -0.5f, -0.5f), new Vector3(0.5f, -0.5f, -0.5f), new Vector3(0.5f, -0.5f, 0.5f), new Vector3(-0.5f, -0.5f, 0.5f) }
    };

    private void GenerateMesh()
    {
        List<Vector3> pos = new List<Vector3>();
        List<Vector3> nor = new List<Vector3>();
        List<Color> col = new List<Color>();
        List<int> idx = new List<int>();

        for (int z = 0; z < zSize; ++z)
        {
            for (int y = 0; y < ySize; ++y)
            {
                for (int x = 0; x < xSize; ++x)
                {
                    Blocks block = (Blocks)blocks[PosToIndex(x, y, z)];

                    if (block == Blocks.AIR)
                    {
                        continue;
                    }

                    var color = block switch
                    {
                        Blocks.STICK => new Color(0.7f, 0.7f, 0.7f),
                        Blocks.CANDY => Color.magenta,
                        _ => Color.grey,
                    };
                    color = OffsetColorRandomly(color);

                    for (int i = 0; i < 6; ++i)
                    {
                        var direction = directions[i];

                        Vector3Int newPos = new Vector3Int(x, y, z) + direction;
                        if (newPos.x >= 0 && newPos.x < xSize && newPos.y >= 0 && newPos.y < ySize && newPos.z >= 0 && newPos.z < zSize)
                        {
                            Blocks neighborBlock = (Blocks)blocks[PosToIndex(newPos.x, newPos.y, newPos.z)];
                            if (neighborBlock != Blocks.AIR)
                            {
                                continue;
                            }
                        }

                        int idx1 = pos.Count;

                        foreach (var vert in directionVerts[i])
                        {
                            pos.Add(vert + new Vector3(x, y, z));
                        }

                        for (int j = 0; j < 4; ++j)
                        {
                            nor.Add(direction);
                            col.Add(color);
                        }

                        idx.Add(idx1);
                        idx.Add(idx1 + 1);
                        idx.Add(idx1 + 2);
                        idx.Add(idx1);
                        idx.Add(idx1 + 2);
                        idx.Add(idx1 + 3);
                    }
                }
            }
        }

        if (voxelMesh != null)
        {
            Destroy(voxelMesh);
        }

        voxelMesh = new Mesh
        {
            name = "Voxel Mesh",
            indexFormat = UnityEngine.Rendering.IndexFormat.UInt32,
            vertices = pos.ToArray(),
            normals = nor.ToArray(),
            colors = col.ToArray(),
            triangles = idx.ToArray()
        };

        GetComponent<MeshFilter>().mesh = voxelMesh;
        if (useCustomShader)
        {
            GetComponent<MeshRenderer>().material = new Material(meshShader);
        }
    }

    private static Color OffsetColorRandomly(Color original, float maxOffset = 0.06f)
    {
        float OffsetValue(float value)
        {
            float offset = Random.Range(-maxOffset, maxOffset);
            return Mathf.Clamp01(value + offset);
        }

        original.r = OffsetValue(original.r);
        original.g = OffsetValue(original.g);
        original.b = OffsetValue(original.b);

        return original;
    }

    private void OnDestroy()
    {
        if (voxelBuffer != null)
        {
            voxelBuffer.Release();
        }

        if (voxelMesh != null)
        {
            Destroy(voxelMesh);
        }
    }
}
